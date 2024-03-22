#include "MFGVRP_Solver.h"

//FRACTIONALITY CUTS::
//USE JENS CUTS TO INCREASE LOWER BOUND



ILOUSERCUTCALLBACK1(FracCuts, MFGVRP_Solver&, solver) {
	//if (getMIPRelativeGap() < 0.015) return;
	double start = 0;
	double end = 0;
	double startTot = solver.Cplex.getTime();
	double StartAddingCuts = 0;
	double endTot = 0;
	double gap = getMIPRelativeGap();
	IloInt i;
	IloEnv env = solver.getEnv();
	IloInt n = solver.getNumCust() + 1;
	int numChargers = solver.getChargers();
	bool printCutInfo = false;
	bool edgesNeedFix = false;
	double edgeVal;
	double edgeValX;
	double edgeValY;
	bool isRoot = false;

	IloExpr expr(env);
	IloExpr expr2(env);
	int* d = new int[n];
	int* s = new int[n];
	int Q = solver.getCapacity();
	int T = solver.getMaxTime();

	std::vector<int> CutList;
	std::vector<int> ExtList;
	bool CutsAdded = false;
	int NNodes = getNnodes();
	int NoOfChargingEdges = 0;
	if (NNodes == 0) isRoot = true;

	int LHS = 0;

	MFGVRP_Solver::UsedEdges EdgesUsed(n);

	std::vector<int> EdgeTail; EdgeTail.push_back(0);
	std::vector<int> EdgeHead; EdgeHead.push_back(0);
	std::vector<double> EdgeX; EdgeX.push_back(-1);
	std::vector<int> unique_ids; unique_ids.push_back(-1);


	int OriginID = 1;
	int Qmin = 0;
	double TotDemand = 0;


	//CVRPSEP parameters:
	CnstrMgrPointer MyCutsCMP;

	CMGR_CreateCMgr(&MyCutsCMP, 100);




	IloNumArray xTemp1d(env, solver.xDummy.getSize());
	IloNumArray yTemp1d(env, solver.yDummy.getSize());

	start = solver.Cplex.getTime();
	getValues(xTemp1d, solver.xDummy);
	getValues(yTemp1d, solver.yDummy);
	end = solver.Cplex.getTime();
	solver.Stats->TotalLoadSepData += end - start;
	start = solver.Cplex.getTime();
	int cntX = 0;
	int cntY = 0;
	double xVal; double yVal;
	for (int i = 0; i < n; i++)
	{
		if (i > 0) { d[i] = solver.getDemand(i); TotDemand += d[i];; s[i] = solver.getServiceTime(i); }


		for (int j = 0; j < n; j++)
		{
			xVal = 0;
			yVal = 0;
			for (int k = 0; k < 2; k++)
			{
				if (xTemp1d[cntX] > 0.001)
				{
					xVal += xTemp1d[cntX];
					if (k == 0)
					{
						EdgesUsed.DirectConnections[i].succ.push_back(j);
						EdgesUsed.DirectConnections[i].succVals.push_back(xTemp1d[cntX]);
						EdgesUsed.DirectConnections[j].pre.push_back(i);
						EdgesUsed.DirectConnections[j].preVals.push_back(xTemp1d[cntX]);
					}
					else
					{
						EdgesUsed.ICEVConnections[i].succ.push_back(j);
						EdgesUsed.ICEVConnections[i].succVals.push_back(xTemp1d[cntX]);
						EdgesUsed.ICEVConnections[j].pre.push_back(i);
						EdgesUsed.ICEVConnections[j].preVals.push_back(xTemp1d[cntX]);
					}
				}

				cntX++;

			}
			for (int r = 0; r < solver.R[i][j].size(); r++)
			{
				if (yTemp1d[cntY] > 0.001)
				{
					xVal += yTemp1d[cntY]; yVal += yTemp1d[cntY];
					EdgesUsed.ChargerConnections[i].succ.push_back(j);
					EdgesUsed.ChargerConnections[i].succVals.push_back(yTemp1d[cntY]);
					EdgesUsed.ChargerConnections[i].succCharger.push_back(solver.R[i][j][r]);
					EdgesUsed.ChargerConnections[j].pre.push_back(i);
					EdgesUsed.ChargerConnections[j].preVals.push_back(yTemp1d[cntY]);
					EdgesUsed.ChargerConnections[j].preCharger.push_back(solver.R[i][j][r]);
				}

				cntY++;
			}
			if (xVal > 0.001) {
				EdgesUsed.Consolidated[i].succ.push_back(j);
				EdgesUsed.Consolidated[i].succVals.push_back(xVal-yVal);
				EdgesUsed.Consolidated[j].pre.push_back(i);
				EdgesUsed.Consolidated[j].preVals.push_back(xVal-yVal);
				if (i > j) {

					auto it = std::find(unique_ids.begin(), unique_ids.end(), 1000 * i + j);
					if (it != unique_ids.end())
					{
						size_t index = std::distance(unique_ids.begin(), it);

						EdgeX[index] += xVal;
					}
					else
					{
						unique_ids.push_back(i * 1000 + j);
						EdgeTail.push_back(j == 0 ? n : j);
						EdgeHead.push_back(i);
						EdgeX.push_back(xVal);

					}

				}
				else {

					auto it = std::find(unique_ids.begin(), unique_ids.end(), 1000 * j + i);
					if (it != unique_ids.end())
					{

						size_t index = std::distance(unique_ids.begin(), it);

						EdgeX[index] += xVal;

					}
					else
					{
						unique_ids.push_back(j * 1000 + i);
						EdgeTail.push_back(i == 0 ? n : i);
						EdgeHead.push_back(j);
						EdgeX.push_back(xVal);
					}

				}
			}
		}
	}



	xTemp1d.end();
	yTemp1d.end();


	end = solver.Cplex.getTime();
	solver.Stats->ConvertVals += end - start;
	start = solver.Cplex.getTime();

	end = solver.Cplex.getTime();
	solver.Stats->SetupSepVectors += end - start;
	Qmin = TotDemand - (n - 1) * Q;
	if (!edgesNeedFix) //If no edges needs to be fixed, we can initiate our seperation problem:
	{
		start = solver.Cplex.getTime();
		int cntCut = MyCutsCMP->Size;


		CAPSEP_SeparateCapCuts(n - 1,
			d,
			Q,
			EdgeTail.size() - 1,
			EdgeTail.data(),
			EdgeHead.data(),
			EdgeX.data(),
			solver.MyOldCutsCMP,
			solver.MaxNoOfCuts,
			solver.EpsForIntegrality,
			&solver.IntegerAndFeasible,
			&solver.MaxViolation,
			MyCutsCMP);



		//FCISEP_SeparateFCIs(n - 1,
		//	d,
		//	Q,
		//	solver.NoOfEdges,
		//	EdgeTail.data(),
		//	EdgeHead.data(),
		//	EdgeX.data(),
		//	solver.MyOldCutsCMP,
		//	20000,
		//	10,
		//	&solver.MaxViolation,
		//	MyCutsCMP);

		if (NNodes == -1 && MyCutsCMP->Size == 0) {

			if (getNiterations() % 2 == 0)
			{
				//if (solver.MaxViolation < 0.05)
				//{
				MSTARSEP_SeparateMultiStarCuts(n - 1,
					d,
					Q,
					EdgeTail.size() - 1,
					EdgeTail.data(),
					EdgeHead.data(),
					EdgeX.data(),
					solver.MyOldCutsCMP,
					solver.MaxNoOfCuts * 2 - MyCutsCMP->Size,
					&solver.MaxViolation,
					MyCutsCMP);
				//}
				if (MyCutsCMP->Size == 0)//&& solver.MaxViolation < 0.1)
				{
					COMBSEP_SeparateCombs(
						n - 1,
						d,
						Q,
						Qmin,
						EdgeTail.size() - 1,
						EdgeTail.data(),
						EdgeHead.data(),
						EdgeX.data(),
						10,
						&solver.MaxViolation,
						MyCutsCMP
					);
				}
			}
			else
			{
				//if ( solver.MaxViolation < 0.1)
				//{
				COMBSEP_SeparateCombs(
					n - 1,
					d,
					Q,
					Qmin,
					EdgeTail.size() - 1,
					EdgeTail.data(),
					EdgeHead.data(),
					EdgeX.data(),
					10,
					&solver.MaxViolation,
					MyCutsCMP
				);
				//}
				if (MyCutsCMP->Size == 0) //&& solver.MaxViolation < 0.05)
				{
					MSTARSEP_SeparateMultiStarCuts(n - 1,
						d,
						Q,
						EdgeTail.size() - 1,
						EdgeTail.data(),
						EdgeHead.data(),
						EdgeX.data(),
						solver.MyOldCutsCMP,
						solver.MaxNoOfCuts * 2 - MyCutsCMP->Size,
						&solver.MaxViolation,
						MyCutsCMP);
				}
			}

		}

		end = solver.Cplex.getTime();
		solver.Stats->CapSepTime += end - start;

		//if (NNodes == 0 && MyCutsCMP->Size == 0)
		//{
		//	solver.ConnectedComponents(&EdgesUsed, MyCutsCMP);
		//}

		//if (MyCutsCMP->Size == 0) {
		//	CAPSEP_SeparateCapCuts(n - 1,
		//		s,
		//		T,
		//		EdgeTail.size() - 1,
		//		EdgeTail.data(),
		//		EdgeHead.data(),
		//		EdgeX.data(),
		//		solver.MyOldCutsCMP,
		//		solver.MaxNoOfCuts,
		//		solver.EpsForIntegrality,
		//		&solver.IntegerAndFeasible,
		//		&solver.MaxViolation,
		//		MyCutsCMP);
		//}



		if (MyCutsCMP->Size == 0) {
			//if (NNodes>0)
			//{
			start = solver.Cplex.getTime();
			if (getNiterations()==417)
			{
				printf("");
			}
			solver.SeparateCuts(&EdgesUsed, MyCutsCMP);
			

			solver.CutCounter = 0;

			if (MyCutsCMP->Size == 0)
			{
				//printf("\nStart connected components");
				solver.ConnectedComponents(&EdgesUsed, MyCutsCMP);
				//printf("\nFinished with %d cuts", MyCutsCMP->Size);
			}



			solver.NodeID = getNodeId();

			end = solver.Cplex.getTime();
			solver.Stats->EnumerationTime += end - start;
		}
		if (MyCutsCMP->Size != 0)
		{
			StartAddingCuts = solver.Cplex.getTime();
			CutsAdded = true;
			for (IloInt cons = 0; cons < MyCutsCMP->Size; cons++)
			{
				CutList.clear();
				ExtList.clear();
				if (MyCutsCMP->CPL[cons]->CType == CMGR_CT_CAP)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->CapCuts++;
						solver.Stats->numberOfCutsAtRootNode->CapCuts++;
					}
					else solver.Stats->totalNumberOfCuts->CapCuts++;

					solver.CAPExpr(expr, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner){
						solver.consDesc.push_back("frac cap");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				//else if (MyCutsCMP->CPL[cons]->CType == CMGR_CT_FCI)
				//{
				//	if (printCutInfo)
				//		printf("\n%d", MyCutsCMP->CPL[cons]->CType);
				//	if (NNodes == 0) {
				//		solver.Stats->totalNumberOfCuts->FCI++;
				//		solver.Stats->numberOfCutsAtRootNode->FCI++;
				//	}
				//	else solver.Stats->totalNumberOfCuts->FCI++;

				//	solver.FCIExpr(expr, MyCutsCMP->CPL[cons]);
				//	add(expr >= MyCutsCMP->CPL[cons]->RHS);
				//}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_CT_STR_COMB)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->Comb++;
						solver.Stats->numberOfCutsAtRootNode->Comb++;
					}
					else solver.Stats->totalNumberOfCuts->Comb++;

					solver.COMBExpr(expr, MyCutsCMP->CPL[cons]);

					add(expr >= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("frac comb");
						solver.cons.add(expr >= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_CT_MSTAR) {
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->MultiStar++;
						solver.Stats->numberOfCutsAtRootNode->MultiStar++;
					}
					else solver.Stats->totalNumberOfCuts->MultiStar++;

					solver.MTSTARExpr(expr, MyCutsCMP->CPL[cons]);
					add(expr >= MyCutsCMP->CPL[cons]->L);
					if (solver.runRefiner) {
						solver.consDesc.push_back("frac MTSTAR");
						solver.cons.add(expr >= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_NO_CHARGE_PATH)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.NoChargePathExpr(expr, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("frac no charge path");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_NO_CHARGE_TOUR)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.NoChargeTourExpr(expr, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("frac no charge tour");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_NO_CHARGE_SET)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.NoChargeSetExpr(expr, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("frac no charge set");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_NO_CHARGE_IN_TOUR)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.ChargerTourExpr(expr, 0, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("frac no charge in tour");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_NO_CHARGE_OUT_TOUR)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.ChargerTourExpr(expr, 1, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("frac no charge out tour");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_NO_CHARGE_DUAL_TOUR)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.ChargerTourExpr(expr, 2, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("frac no charge dual tour");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_FIXED_IN_CHARGER)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.FixedChargerExpr(expr, 0, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("frac fixed in charger");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_FIXED_OUT_CHARGER)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.FixedChargerExpr(expr, 1, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("frac fixed out charger");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_FIXED_CHARGERS)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.FixedChargerExpr(expr, 2, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("frac fixed chargers");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_FIXED_IN_EDGE)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					if (MyCutsCMP->CPL[cons]->IntListSize==5)
					{
						if (MyCutsCMP->CPL[cons]->IntList[2] == 25 && MyCutsCMP->CPL[cons]->IntList[3] == 9 && MyCutsCMP->CPL[cons]->IntList[4] == 1 && MyCutsCMP->CPL[cons]->IntList[5] == 12)
							printf("Hey");
					}
					
					solver.FixedChargingEdgeExpr(expr, 0, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("frac fixed in charging edge");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_FIXED_OUT_EDGE)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.FixedChargingEdgeExpr(expr, 1, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("frac fixed out charging edge");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_FIXED_EDGES)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.FixedChargingEdgeExpr(expr, 2, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("frac fixed charging edges");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_DURATION_TOUR_BOTH)
				{
				if (printCutInfo)
					printf("\n%d", MyCutsCMP->CPL[cons]->CType);
				solver.DurationTourExpr(expr, 2, MyCutsCMP->CPL[cons]);
				add(expr <= MyCutsCMP->CPL[cons]->RHS);
				if (solver.runRefiner) {
					solver.consDesc.push_back("frac duration tour both");
					solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
				}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_DURATION_TOUR_EV)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.DurationTourExpr(expr, 0, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("frac duration tour EV");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_DURATION_TOUR_ICEV)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.DurationTourExpr(expr, 1, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("frac duration tour ICEV");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_DURATION_SET_BOTH)
				{
				if (printCutInfo)
					printf("\n%d", MyCutsCMP->CPL[cons]->CType);
				solver.DurationSetExpr(expr, 2, MyCutsCMP->CPL[cons]);
				add(expr <= MyCutsCMP->CPL[cons]->RHS);
				if (solver.runRefiner) {
					solver.consDesc.push_back("frac duration set");
					solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
				}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_DURATION_SET_EV)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.DurationSetExpr(expr, 0, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("frac duration set EV");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_DURATION_SET_ICEV)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.DurationSetExpr(expr, 1, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("frac duration set ICEV");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
			}


			for (i = 0; i < MyCutsCMP->Size; i++)
			{
				CMGR_MoveCnstr(MyCutsCMP, solver.MyOldCutsCMP, i, 0);
			}

			MyCutsCMP->Size = 0;
			end = solver.Cplex.getTime();
			solver.Stats->AddCutsTime += end - StartAddingCuts;

		}
	}
	delete[] d;
	delete[] s;

	if (NNodes == 0) {
		solver.Stats->LBAtRootNode = getBestObjValue();
		solver.Stats->timeRootNode = solver.Cplex.getTime();
	}

	if (solver.Stats->treeDepth < getCurrentNodeDepth()) solver.Stats->treeDepth = getCurrentNodeDepth();
	endTot = solver.Cplex.getTime();
	solver.Stats->TotalSepTime += endTot - startTot;
}
// LAZY CUTS::
//RIGHT NOW CAP CUTS ADDED TOGETHER WITH BATTERY CUTS. CAP CUTS WORKS AS INTENDED, BATTERY CUTS NEEDS TO BE CHECKED FOR WHETHER THEY ALWAYS REACH CHARGING STATION BEFORE RUNNING OUT OF BATTERY:
ILOLAZYCONSTRAINTCALLBACK1(LazyCapCuts, MFGVRP_Solver&, solver) {
	double start = 0;
	double end = 0;
	double startTot = solver.Cplex.getTime();
	double endTot = 0;
	bool printCutInfo = false;
	double startAddingCuts = 0;
	IloInt i;
	IloEnv env = solver.getEnv();
	IloInt n = solver.x.getSize();
	bool edgesNeedFix = false;


	double edgeVal;
	double edgeValX;
	double edgeValY;
	IloExpr expr(env);

	int* d = new int[n];
	int Q = solver.getCapacity();
	int* s = new int[n];
	int T = solver.getMaxTime();

	bool CutsAdded = false;

	std::vector<int> CutList;
	std::vector<int> ExtList;
	int NNodes = getNnodes();


	int LHS = 0;

	int NoOfChargingEdges = 0;

	MFGVRP_Solver::UsedEdges EdgesUsed(n);

	std::vector<int> EdgeTail; EdgeTail.push_back(0);
	std::vector<int> EdgeHead; EdgeHead.push_back(0);
	std::vector<double> EdgeX; EdgeX.push_back(-1);
	std::vector<int> unique_ids; unique_ids.push_back(-1);

	//CVRPSEP parameters:
	CnstrMgrPointer MyCutsCMP;

	CMGR_CreateCMgr(&MyCutsCMP, 100);


	IloNumArray xTemp1d(env, solver.xDummy.getSize());
	IloNumArray yTemp1d(env, solver.yDummy.getSize());

	start = solver.Cplex.getTime();
	getValues(xTemp1d, solver.xDummy);
	getValues(yTemp1d, solver.yDummy);
	end = solver.Cplex.getTime();
	solver.Stats->TotalLoadSepData += end - start;
	start = solver.Cplex.getTime();

	int Qmin = 0;
	int cntX = 0;
	int cntY = 0;
	double xVal; double yVal;


	for (int i = 0; i < n; i++)
	{
		if (i > 0) { d[i] = solver.getDemand(i); Qmin += d[i];; s[i] = solver.getServiceTime(i); }

		for (int j = 0; j < n; j++)
		{
			xVal = 0;
			yVal = 0;
			//xTemp[i][j] = IloNumArray(env, 2);
			for (int k = 0; k < 2; k++)
			{
				//xTemp[i][j][k] = xTemp1d[cntX];	
				if (xTemp1d[cntX] > 0.001)
				{
					xVal += xTemp1d[cntX];
					if (k == 0)
					{
						EdgesUsed.DirectConnections[i].succ.push_back(j);
						EdgesUsed.DirectConnections[i].succVals.push_back(xTemp1d[cntX]);
						EdgesUsed.DirectConnections[j].pre.push_back(i);
						EdgesUsed.DirectConnections[j].preVals.push_back(xTemp1d[cntX]);
					}
					else
					{
						EdgesUsed.ICEVConnections[i].succ.push_back(j);
						EdgesUsed.ICEVConnections[i].succVals.push_back(xTemp1d[cntX]);
						EdgesUsed.ICEVConnections[j].pre.push_back(i);
						EdgesUsed.ICEVConnections[j].preVals.push_back(xTemp1d[cntX]);
					}
				}

				cntX++;

			}
			for (int r = 0; r < solver.R[i][j].size(); r++)
			{
				//yTemp[i][j][r] = yTemp1d[cntY];
				if (yTemp1d[cntY] > 0.001)
				{
					xVal += yTemp1d[cntY]; yVal += yTemp1d[cntY];
					EdgesUsed.ChargerConnections[i].succ.push_back(j);
					EdgesUsed.ChargerConnections[i].succVals.push_back(yTemp1d[cntY]);
					EdgesUsed.ChargerConnections[i].succCharger.push_back(solver.R[i][j][r]);
					EdgesUsed.ChargerConnections[j].pre.push_back(i);
					EdgesUsed.ChargerConnections[j].preVals.push_back(yTemp1d[cntY]);
					EdgesUsed.ChargerConnections[j].preCharger.push_back(solver.R[i][j][r]);
				}

				cntY++;
			}
			if (xVal > 0.001) {
				EdgesUsed.Consolidated[i].succ.push_back(j);
				EdgesUsed.Consolidated[i].succVals.push_back(xVal-yVal);
				EdgesUsed.Consolidated[j].pre.push_back(i);
				EdgesUsed.Consolidated[j].preVals.push_back(xVal-yVal);
				if (i > j) {
					auto it = std::find(unique_ids.begin(), unique_ids.end(), 1000 * i + j);
					if (it != unique_ids.end())
					{
						size_t index = std::distance(unique_ids.begin(), it);

						EdgeX[index] += xVal;
					}
					else
					{
						unique_ids.push_back(i * 1000 + j);
						EdgeTail.push_back(j == 0 ? n : j);
						EdgeHead.push_back(i);
						EdgeX.push_back(xVal);

					}

				}
				else {

					auto it = std::find(unique_ids.begin(), unique_ids.end(), 1000 * j + i);
					if (it != unique_ids.end())
					{

						size_t index = std::distance(unique_ids.begin(), it);

						EdgeX[index] += xVal;
					}
					else
					{
						unique_ids.push_back(j * 1000 + i);
						EdgeTail.push_back(i == 0 ? n : i);
						EdgeHead.push_back(j);
						EdgeX.push_back(xVal);
					}

				}
			}
		}
	}


	xTemp1d.end();
	yTemp1d.end();
	end = solver.Cplex.getTime();


	if (!edgesNeedFix) //If no edges needs to be fixed, we can initiate our seperation problem:
	{
		start = solver.Cplex.getTime();
		int cntCut = MyCutsCMP->Size;
		CAPSEP_SeparateCapCuts(n - 1,
			d,
			Q,
			EdgeTail.size() - 1,
			EdgeTail.data(),
			EdgeHead.data(),
			EdgeX.data(),
			solver.MyOldCutsCMP,
			solver.MaxNoOfCuts,
			solver.EpsForIntegrality,
			&solver.IntegerAndFeasible,
			&solver.MaxViolation,
			MyCutsCMP);

		//if (MyCutsCMP->Size == 0)
		//{
		//	CAPSEP_SeparateCapCuts(n - 1,
		//		s,
		//		T,
		//		EdgeTail.size() - 1,
		//		EdgeTail.data(),
		//		EdgeHead.data(),
		//		EdgeX.data(),
		//		solver.MyOldCutsCMP,
		//		solver.MaxNoOfCuts,
		//		solver.EpsForIntegrality,
		//		&solver.IntegerAndFeasible,
		//		&solver.MaxViolation,
		//		MyCutsCMP);
		//}
		end = solver.Cplex.getTime();
		solver.Stats->CapSepTime += end - start;


		if (MyCutsCMP->Size == 0) {

			start = solver.Cplex.getTime();
			solver.SeparateCuts(&EdgesUsed, MyCutsCMP);
			solver.CutCounter = 0;

			if (MyCutsCMP->Size == 0)
			{
				//printf("\nStart connected components");
				solver.ConnectedComponents(&EdgesUsed, MyCutsCMP);
				//printf("\nFinished with %d cuts", MyCutsCMP->Size);
			}

			if (MyCutsCMP->Size == 0)
			{
				solver.FinalFeasibilityCheck(&EdgesUsed.DirectConnections, &EdgesUsed.ChargerConnections, &EdgesUsed.ICEVConnections, MyCutsCMP);
			}

			//if (MyCutsCMP->Size == 0 && NNodes==0)
			//{
			//	printf("\nStart connected components");
			//	solver.ConnectedComponents(&EdgesUsed, MyCutsCMP);
			//	printf("\nFinished with %d cuts", MyCutsCMP->Size);
			//}

			end = solver.Cplex.getTime();
			solver.Stats->EnumerationTime += end - start;
		}//, EnergyCuts);




		//std::cout << "\n\nNumber of violated cuts: " << MyCutsCMP->Size << "\n\nCuts are:";
		if (MyCutsCMP->Size != 0)
		{
			startAddingCuts = solver.Cplex.getTime();
			CutsAdded = true;
			for (IloInt cons = 0; cons < MyCutsCMP->Size; cons++)
			{
				CutList.clear();
				ExtList.clear();
				if (MyCutsCMP->CPL[cons]->CType == CMGR_CT_CAP)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->CapCuts++;
						solver.Stats->numberOfCutsAtRootNode->CapCuts++;
					}
					else solver.Stats->totalNumberOfCuts->CapCuts++;

					solver.CAPExpr(expr, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("int cap");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				//else if (MyCutsCMP->CPL[cons]->CType == CMGR_CT_FCI)
				//{
				//	if (printCutInfo)
				//		printf("\n%d", MyCutsCMP->CPL[cons]->CType);
				//	if (NNodes == 0) {
				//		solver.Stats->totalNumberOfCuts->FCI++;
				//		solver.Stats->numberOfCutsAtRootNode->FCI++;
				//	}
				//	else solver.Stats->totalNumberOfCuts->FCI++;

				//	solver.FCIExpr(expr, MyCutsCMP->CPL[cons]);
				//	add(expr >= MyCutsCMP->CPL[cons]->RHS);
				//}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_CT_STR_COMB)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->Comb++;
						solver.Stats->numberOfCutsAtRootNode->Comb++;
					}
					else solver.Stats->totalNumberOfCuts->Comb++;

					solver.COMBExpr(expr, MyCutsCMP->CPL[cons]);

					add(expr >= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("int comb");
						solver.cons.add(expr >= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_CT_MSTAR) {
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->MultiStar++;
						solver.Stats->numberOfCutsAtRootNode->MultiStar++;
					}
					else solver.Stats->totalNumberOfCuts->MultiStar++;

					solver.MTSTARExpr(expr, MyCutsCMP->CPL[cons]);
					add(expr >= MyCutsCMP->CPL[cons]->L);
					if (solver.runRefiner) {
						solver.consDesc.push_back("int MTSTAR");
						solver.cons.add(expr >= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_NO_CHARGE_PATH)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.NoChargePathExpr(expr, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("frac no charge path");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_NO_CHARGE_TOUR)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.NoChargeTourExpr(expr, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("int no charge tour");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_NO_CHARGE_SET)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.NoChargeSetExpr(expr, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("int no charge set");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_NO_CHARGE_IN_TOUR)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.ChargerTourExpr(expr, 0, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("int no charge in tour");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_NO_CHARGE_OUT_TOUR)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.ChargerTourExpr(expr, 1, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("int no charge out tour");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_NO_CHARGE_DUAL_TOUR)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.ChargerTourExpr(expr, 2, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("int no charge dual tour");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_FIXED_IN_CHARGER)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.FixedChargerExpr(expr, 0, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("int fixed in charger");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_FIXED_OUT_CHARGER)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.FixedChargerExpr(expr, 1, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("int fixed out charger");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_FIXED_CHARGERS)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.FixedChargerExpr(expr, 2, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("int fixed chargers");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_FIXED_IN_EDGE)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.FixedChargingEdgeExpr(expr, 0, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("int fixed in charging edge");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_FIXED_OUT_EDGE)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.FixedChargingEdgeExpr(expr, 1, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("int fixed out charging edge");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_FIXED_EDGES)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.FixedChargingEdgeExpr(expr, 2, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("int fixed charging edges");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_DURATION_TOUR_BOTH)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.DurationTourExpr(expr, 2, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("int duration tour both");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_DURATION_TOUR_EV)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.DurationTourExpr(expr, 0, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("int duration tour EV");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_DURATION_TOUR_ICEV)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.DurationTourExpr(expr, 1, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("int duration tour ICEV");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_DURATION_SET_BOTH)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.DurationSetExpr(expr, 2, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("int duration set");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_DURATION_SET_EV)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.DurationSetExpr(expr, 0, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("int duration set EV");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_DURATION_SET_ICEV)
				{
					if (printCutInfo)
						printf("\n%d", MyCutsCMP->CPL[cons]->CType);
					solver.DurationSetExpr(expr, 1, MyCutsCMP->CPL[cons]);
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					if (solver.runRefiner) {
						solver.consDesc.push_back("int duration set ICEV");
						solver.cons.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					}
				}
			}

		}
	}

	delete[] d;
	delete[] s;



	if (getNnodes() == 0) {
		solver.Stats->LBAtRootNode = getBestObjValue();
		solver.Stats->timeRootNode = solver.Cplex.getTime();
	}

	if (solver.Stats->treeDepth < getCurrentNodeDepth()) solver.Stats->treeDepth = getCurrentNodeDepth();

	if (!CutsAdded) {
		if (solver.Stats->treeFeasDepth == -1) solver.Stats->treeFeasDepth = getCurrentNodeDepth();
		if (solver.Stats->currentUB > getObjValue()) {
			solver.Stats->currentUB = getObjValue();
			solver.Stats->treeUBDepth = getCurrentNodeDepth();
		}
	}
	endTot = solver.Cplex.getTime();
	solver.Stats->TotalSepTime += endTot - startTot;
}
//TSP cuts LAZY:
ILOLAZYCONSTRAINTCALLBACK1(LazyTSPCuts, MFGVRP_Solver&, solver) {
	
	try
	{

	
		IloEnv env = solver.getEnv();
		IloExpr expr(env);
	


		double ObjVal = getObjValue();
		double B = solver.getBatteryCap();
		double T = solver.getMaxTime();

		int NoOFCustomers = solver.getNoTSPCustomers();
		int* TSPCustomers = new int[NoOFCustomers];
		for (int i = 0; i < NoOFCustomers; i++) TSPCustomers[i] = solver.getTSPCutomers(i);
		int idx = TSPCustomers[0];
		int cnt = 0;
		std::vector<int> SubTour;
		if (getBestObjValue()<=solver.TSPCutOff+0.03)
		{
			int f = solver.getChargers();
			int n = solver.getNumCust();
			IloNum2DMatrix xVals(env, n + f + 1);


			for (IloInt i = 0; i < n + f + 1; i++) {
				xVals[i] = IloNumArray(env, n + f + 1);
				getValues(xVals[i], solver.xTSP[i]);
			}
			while (true)
			{
				for (int i = 0; i < NoOFCustomers; i++)
				{
					if (xVals[idx][TSPCustomers[i]] >= 0.99) {

						idx = TSPCustomers[i];
						SubTour.push_back(idx);
						if (idx != TSPCustomers[0]) { ++cnt; }
						break;
					}
				}
				if (idx == TSPCustomers[0]) {
					if (cnt < NoOFCustomers - 1)
					{
						expr.clear();
						for (int i = 0; i < SubTour.size(); i++) for (int j = 0; j < SubTour.size(); j++) expr += solver.xTSP[SubTour[i]][SubTour[j]];
						add(expr <= cnt);
						break;
					}
					else break;
				}
			}
			for (int i = 0; i < n + f + 1; i++)
				xVals[i].end();
			xVals.end();
		}


		delete[] TSPCustomers;
	}
	catch (const std::exception&)
	{
		printf("\nSOMETHING FISHY IS GOING ON!!!!!!\n");
	}
}
//Define constructor:
MFGVRP_Solver::MFGVRP_Solver() {


	CMGR_CreateCMgr(&MyOldCutsCMP, 100);

	//Initializing CPLEX models:
	Model = IloModel(env);

	Cplex = IloCplex(Model);
	//Cplex.setParam(IloCplex::WorkDir, "c:/cplex/");
	//Cplex.setParam(IloCplex::NodeFileInd, 3);
	//Cplex.setParam(IloCplex::Param::MIP::Tolerances::UpperCutoff, 309.5);
	//Cplex.setParam(IloCplex::WorkMem, 500.0);
	//Cplex.setParam(IloCplex::SolnPoolCapacity, 1);
	//Cplex.setParam(IloCplex::Param::MIP::Limits::Nodes, 0); //Only process root node
	//Cplex.setParam(IloCplex::Threads, 1);
	//Cplex.setParam(IloCplex::ParallelMode, 1);
	//Cplex.setParam(IloCplex::EpRHS, 1e-9);
	//Cplex.setParam(IloCplex::EpAGap, 0.0);
	Cplex.setParam(IloCplex::VarSel, 4);
	Cplex.setParam(IloCplex::MemoryEmphasis, true);
	//Variables
	x = IloVar3DMatrix(env);
	y = IloVar3DMatrix(env);
	//SOC = IloNumVarArray(env);
	//e = IloVar2DMatrix(env);
	//time = IloNumVarArray(env);

	//Parameters:
	c = IloNum3DMatrix(env); //Cost for driving between customer i and j with charging mode p
	ICEV_Cost_matrix = IloNum2DMatrix(env); //Cost for driving between customer i and j with charging mode p
	d = IloNumArray(env); //Customer demands
	s = IloNumArray(env); //Contains service times all customers
	DrawCounter = 0;

	lazy = IloRangeArray(env);
	user = IloRangeArray(env);
	//infeas = IloRangeArray(env);
	//Constraints:
	cons = IloRangeArray(env); //Holds all initial constraints
	//infeas = IloRangeArray(env);
	CapCuts = IloRangeArray(env); //Holds all capacity cuts
	//Expressions
	expr = IloExpr(env); //New expression 1
	expr2 = IloExpr(env); //New expression 2
	//infeas = IloRangeArray(env);

	//Fractional cuts:
	EpsForIntegrality = 0.0001;
	MaxViolation = 0;

	//Statistics:
	Stats = new testStats;
}

//---------------------------------------------------------------------------------------------------------------------------//
MFGVRP_Solver::~MFGVRP_Solver() {
	if (!HasCleanedUp) {
		//Freeing memory occupied by variables:
		if (x.getImpl()) {
			if (n > 0) for (IloInt i = 0; i < x.getSize(); i++) x[i].end();
			x.end();
		}
		if (y.getImpl()) {
			for (int i = 0; i < y.getSize(); i++)
			{
				for (int j = 0; j < y[i].getSize(); j++)
				{
					y[i][j].end();
				}
				y[i].end();
			}
			y.end();
		}

		if (c.getImpl()) {
			for (int i = 0; i < c.getSize(); i++) {
				for (int j = 0; j < c.getSize(); j++)
					c[i][j].end();
				c[i].end();
			}
			c.end();
		}

		if (ICEV_Cost_matrix.getImpl())
		{
			for (int i = 0; i < ICEV_Cost_matrix.getSize(); i++)
				ICEV_Cost_matrix[i].end();
			ICEV_Cost_matrix.end();
		}
		if (d.getImpl()) d.end();
		if (s.getImpl()) s.end();
		if (SOC.getImpl()) SOC.end();
		if (time.getImpl()) time.end();
		if (e.getImpl())
		{
			for (int i = 0; i < e.getSize(); i++)
			{
				e[i].end();
			}
			e.end();

		}
		if (env.getImpl()) env.end();
		//if (Cplex.getImpl()) Cplex.end();
		//if (cons.getImpl()) cons.end();
		//if (Model.getImpl()) Model.end();
		//
		//if (TSPCplex.getImpl()) TSPModel.end();
		//if (TSPCons.getImpl()) TSPCons.end();
		//if (TSPModel.getImpl()) TSPModel.end();



		//Freeing memory occupied by statistics:
		FreeMemoryFromTestStats(Stats);
		delete Stats;
		HasCleanedUp = true;
	}
}

//---------------------------------------------------------------------------------------------------------------------------//
//Define load data
void MFGVRP_Solver::LoadData(const std::string &FileName, bool RunAllCuts) {
	RunCuts = RunAllCuts;
	double SomeNumber; //Used to load in data to cplex arrays
	std::string SomeValue;
	double x_in, y_in;
	std::vector<double> line;
	std::ifstream in(FileName, std::ios::app);
	if (!in) {
		std::cout << "Cannot open file.\n";
	}
	runRefiner = false;
	in >> n; //Number of customers 
	in >> f; //Number of charging stations
	in >> E; //Number of electric vehicles
	in >> C; //Number of conventional vehicles
	in >> G; //Number of green zones
	in >> T; //Max tour duration
	in >> B; //Battery capacity 
	in >> g; //consumption rate
	in >> r; //Refueling rate
	in >> Q; //Load capacity


	//Load in demands:
	for (IloInt i = 0; i < n + 1; i++)
	{
		in >> SomeNumber;
		d.add(SomeNumber);
	}

	//Load in service times:
	for (IloInt i = 0; i < n + 1; i++)
	{
		in >> SomeNumber;
		s.add(SomeNumber);
	}

	//std::cout << "\n\nC_DIR_ij: ";
	//for (IloInt i = 0; i < n + f + 1; i++)
	//{
	//	std::cout << "\n";
	//	for (IloInt j = 0; j < n + 1; j++)
	//	{
	//		in >> SomeNumber;
	//		std::cout << SomeNumber << " ";
	//	}
	//}

	//Load coordinates
	for (IloInt i = 0; i < n + f + 1; i++)
	{
		in >> SomeNumber;
		xCoord.push_back(SomeNumber);
		in >> SomeNumber;
		yCoord.push_back(SomeNumber);
	}

	in>>GZType;
	in >> GZPenalty;

	////Save customer points:
	////for (IloInt i = 1; i < n + 1; i++) CustomersPoints.push_back(std::make_pair(xCoord[i], yCoord[i]));
	//in >> GZType;
	//Load green zone points:
	for (int i = 0; i < G; i++)
	{
		in >> SomeNumber;
		
		for (int j = 0; j < SomeNumber; j++)
		{
			line.clear();
			in >> x_in;
			in >> y_in;
			line.push_back(x_in);
			line.push_back(y_in);
			in >> x_in;
			in >> y_in;
			line.push_back(x_in);
			line.push_back(y_in);
			GreenZones.push_back(line);
		}
	}

	c = IloNum3DMatrix(env, n + f + 1);
	for (IloInt i = 0; i < n + f + 1; i++)
	{
		c[i] = IloNum2DMatrix(env, n + f + 1); //Add dimension
		//std::cout << "\n";
		for (IloInt j = 0; j < n + f + 1; j++)
		{
			c[i][j] = IloNumArray(env, 2);
			in >> c[i][j][0];
			//printf("\n%f", c[i][j][0]);
			//c[i][j] = calc_distances(xCoord[i] - xCoord[j], yCoord[i] - yCoord[j]);
			//std::cout << c[i][j] << " ";
		}
	}

	for (IloInt i = 0; i < n + f + 1; i++)
		for (IloInt j = 0; j < n + f + 1; j++)
			in >> c[i][j][1];


	ICEV_Cost_matrix = IloNum2DMatrix(env, n + 1);
	for (int i = 0; i < n + 1; i++){
		ICEV_Cost_matrix[i] = IloNumArray(env, n + 1); //Add dimension
		for (int j = 0; j < n+1; j++)
		{
			if (GZType>0)
			{
				in >> SomeValue;

				if (SomeValue == "inf")
				{
					ICEV_Cost_matrix[i][j] = 999999;
				}
				else
				{
					ICEV_Cost_matrix[i][j] = std::stod(SomeValue);
				}
			}
			else ICEV_Cost_matrix[i][j] = c[i][j][1];
			
			
		}
	}
	//printf("\n\nEV resource consumption");
	//for (int i = 0; i < c.getSize(); i++)
	//{
	//	printf("\n");
	//	for (int j = 0; j < c[i].getSize(); j++)
	//		printf("%f \t",c[i][j][0]);
	//}

	//printf("\n\nICEV resource consumption");
	//for (int i = 0; i < c.getSize(); i++)
	//{
	//	printf("\n");
	//	for (int j = 0; j < c[i].getSize(); j++)
	//		printf("%f \t", c[i][j][1]);
	//}

	//printf("\n\nICEV cost");
	//for (int i = 0; i < ICEV_Cost_matrix.getSize(); i++)
	//{
	//	printf("\n");
	//	for (int j = 0; j < ICEV_Cost_matrix[i].getSize(); j++)
	//		printf("%f \t", ICEV_Cost_matrix[i][j]);
	//}
	//std::cout << "\nDistance to charging stations: ";
	////Load distances to charging nodes
	//for (IloInt i = 0; i < f; i++)
	//{
	//	u.add(IloNumArray(env, n + 1));
	//	std::cout << "\n";
	//	for (IloInt j = 0; j < n + 1; j++)
	//	{
	//		in >> c[i][j];
	//		std::cout << c[i][j] << " ";
	//	}
	//}
	//std::cout << "\nc_REC_ij: ";
	////Find minimum charging paths
	//for (IloInt i = 0; i < n + 1; i++)
	//{
	//	std::cout << "\n";
	//	for (IloInt j = 0; j < n + 1; j++)
	//	{
	//		minDist = 1000;
	//		for (IloInt p = 0; p < f; p++)
	//		{
	//			if (minDist > u[p][i] + u[p][j])
	//			{
	//				minDist = u[p][i] + u[p][j];
	//				MinDistCharg= p;
	//			}
	//		}
	//		if (i != j) {c[i][j][1] = minDist;}
	//		else { c[i][j][1] = 999;}
	//		std::cout << c[i][j][1] << " ";
	//		
	//	}
	//}

	in.close();
	//SetupGreenZone();
	int numY = NonDominatedChargingPaths();

	//Setup variables:
	y = IloVar3DMatrix(env, n + 1); //Setting 1st dimensions for y variables
	x = IloVar3DMatrix(env, n + 1); //Setting 1st dimensions for x variables
	xDummy = IloNumVarArray(env, (n + 1) * (n + 1) *2, 0, 1, ILOINT);
	yDummy = IloNumVarArray(env, numY, 0, 1, ILOINT);

	for (IloInt i = 0; i < n + 1; i++)
	{
		x[i] = IloVar2DMatrix(env, n + 1); //Setting 2nd dimensions for x variables
		y[i] = IloVar2DMatrix(env, n + 1); //Setting 2nd dimensions for y variables
		for (IloInt j = 0; j < n + 1; j++)
		{
			x[i][j] = IloNumVarArray(env, 2, 0, 1, ILOINT); 
			y[i][j] = IloNumVarArray(env, R[i][j].size(), 0, 1, ILOINT); 
			//xSol1[i][j] = IloNumArray(env, 2);
			for (IloInt k = 0; k < 2; k++)
			{
				std::string xName = "x_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k);
				x[i][j][k].setName(xName.c_str());
			}
			for (IloInt r = 0; r < R[i][j].size(); r++)
			{
				std::string yName = "y_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(r);
				y[i][j][r].setName(yName.c_str());
			}
		}
	}

	//Initiate vectors for identifying paths:
	//AddedToPath = new bool[n + f + 1];
	//LevelNode = new int[n + f + 1];
	//LevelIdx = new int[n + f + 1];
	//Slack = new double[n + f + 1];
	//SlackSet = new double[n + f + 1];
	//Esum = new double[n + f + 1];
	//Tsum = new double[n + f + 1];
	//rhs=IloIntArray (env, n * 2 + 2*f + 2);

	//Set number of cap cuts:
	if(n>50)MaxNoOfCuts = 10; else MaxNoOfCuts = 5;

	Stats->instanceName = FileName;
	Stats->n = n;
	Stats->r = f;

	BuildModels();
}

//---------------------------------------------------------------------------------------------------------------------------//
void MFGVRP_Solver::BuildModels() {

	int count = 0;
	int consindex = 0;
	std::string consName = "";
	//PROBLEM SETUP
	//Set tolerance:
	Cplex.setParam(IloCplex::EpGap, 0);

	//Build objective function: (1)
	expr.clear();
	for (IloInt i = 0; i < n + 1; i++)
	{
		for (IloInt j = 0; j < n + 1; j++)
		{

			expr+= c[i][j][0]* x[i][j][0];
			expr +=ICEV_Cost_matrix[i][j] * x[i][j][1];
			for (IloInt r = 0; r < R[i][j].size(); r++)
			{
				expr+= (c[R[i][j][r]][i][0] + c[R[i][j][r]][j][0])*y[i][j][r];
			}

		}
	}

	Obj = IloObjective(env);
	Obj = IloMinimize(env,expr);
	Model.add(Obj);
	//Cons 2: All customers must be visited by either ICEV or EV (2)
	for (IloInt j = 1; j < n + 1; j++)
	{
		expr.clear();
		//expr2.clear();
		for (IloInt i = 0; i < n + 1; i++)
		{
			for (IloInt k = 0; k < 2; k++)
			{
				expr += x[i][j][k];
				//expr2 += x[j][i][k];
			}
			for (IloInt r = 0; r < R[i][j].size(); r++)
			{
				expr += y[i][j][r];
				//expr2 += y[j][i][r];
			}
		}
		cons.add(expr == 1);
		if(runRefiner)
			consDesc.push_back("degree");
		consName = "Must_enter_customer_" + std::to_string(j);
		cons[count].setName(consName.c_str());
		count++;
	}

	expr.clear();
	expr2.clear();


	for (int i = 0; i < n + 1; i++)
	{
		for (int k = 0; k < 2; k++) x[i][i][k].setUB(0);
		for (int r = 0; r < R[i][i].size(); r++) y[i][i][r].setUB(0);
		for (int r = 0; r < R[0][i].size(); r++) if (R[0][i][r] == n + 1) { y[0][i][r].setUB(0); break; }
		for (int r = 0; r < R[i][0].size(); r++) if (R[i][0][r] == n + 1) { y[i][0][r].setUB(0); break; }
	}


	// Cons 3+4: If a vehicle enter a node, it must leave it again 
	for (IloInt j = 0; j < n + 1; j++)
	{
		expr.clear(); //(3)
		for (IloInt i = 0; i < n + 1; i++)
		{
			expr += x[i][j][0] - x[j][i][0];
			for (IloInt r = 0; r < R[i][j].size(); r++)
			{
				expr += y[i][j][r] - y[j][i][r];
			}
		}
		cons.add(expr == 0);
		if (runRefiner)
			consDesc.push_back("flow EV");
		//infeas.add(expr == 0);
		consName = "Must_leave_customer_" + std::to_string(j) + "_if_visited_by_EV";
		cons[count].setName(consName.c_str());
		count++;
		expr.clear(); //(4)
		for (IloInt i = 0; i < n + 1; i++)
		{
			expr += x[i][j][1] - x[j][i][1];
		}
		cons.add(expr == 0);
		//infeas.add(expr == 0);
		if (runRefiner)
			consDesc.push_back("flow ICEV");
		consName = "Must_leave_customer_" + std::to_string(j) + "_if_visited_by_ICEV";
		cons[count].setName(consName.c_str());
		count++;

	}
	//Setup dummy variables for solution fetching 
	int cnt = 0;
	
	for (int i = 0; i < n + 1; i++)
	{
		for (int j = 0; j < n + 1; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				cons.add(x[i][j][k] - xDummy[cnt] == 0);
				if (runRefiner)
					consDesc.push_back("dummy");
				cnt++;
			}
		}
	}


	cnt = 0;

	for (int i = 0; i < n + 1; i++)
	{
		for (int j = 0; j < n + 1; j++)
		{
			for (int r = 0; r < R[i][j].size(); r++)
			{
				cons.add(y[i][j][r] - yDummy[cnt] == 0);
				if (runRefiner)
					consDesc.push_back("dummy");
				cnt++;
			} 
		}
	}


	Model.add(cons);
	Cplex.setParam(IloCplex::TiLim, 10800);

	Run();

}

//---------------------------------------------------------------------------------------------------------------------------//
void MFGVRP_Solver::Run() {
	Cplex.use(LazyCapCuts(env, *this));
	Cplex.use(FracCuts(env, *this));
	double startTime = Cplex.getTime();

	try
	{
		if (Cplex.solve()) {
			Stats->bestUpperBound = Cplex.getObjValue();
			for (int i = 0; i < n + 1; i++)
			{
				Stats->EVs += Cplex.getValue(x[0][i][0]);
				Stats->ICEVs += Cplex.getValue(x[0][i][1]);
				for (int r = 0; r < R[0][i].size(); r++)
				{
					Stats->EVs += Cplex.getValue(y[0][i][r]);
				}
			}


		}




		
		if (runRefiner && Cplex.getObjValue()> 497.5){


			IloModel FeasModel(env); //TSP model
			IloCplex FeasCplex(FeasModel); //TSP problem solver
			
			FeasModel.add(Obj);
			FeasModel.add(cons);

			FeasCplex.solve();
			std::cout << "\n" << "Cplex stats:\n" << FeasCplex.getStatus() << "\n";
			IloNumArray preferences(env);
			IloConstraintArray Infeasible(env);
			Infeasible.add(cons);
			for (IloInt i = 0; i < Infeasible.getSize(); i++) {
				preferences.add(1.0);  // user may wish to assign unique preferences
			}

			if (FeasCplex.refineConflict(Infeasible, preferences)) {
				IloCplex::ConflictStatusArray conflict = FeasCplex.getConflict(Infeasible);
				env.getImpl()->useDetailedDisplay(IloTrue);
				std::cout << "Conflict :" << std::endl;
				for (IloInt i = 0; i < Infeasible.getSize(); i++) {
					if (conflict[i] == IloCplex::ConflictMember)
						std:: cout << "Proved  : " <<consDesc[i]<< ": " << Infeasible[i] << std::endl;
					if (conflict[i] == IloCplex::ConflictPossibleMember)
						std::cout << "Possible: " << Infeasible[i] << std::endl;
				}
			}
			else
				std::cout << "Conflict could not be refined" << std::endl;
			std::cout << std::endl;
		}

		Stats->time =  Cplex.getTime() - startTime;
		Stats->treeSize = Cplex.getNnodes();

		if (Cplex.getStatus() == IloAlgorithm::Status::Optimal)
			Stats->isOptimal = true;
		else
			Stats->isOptimal = false;
		if (Cplex.getStatus() == IloAlgorithm::Status::Infeasible)
			Stats->isInfeasible = true;
		else
			Stats->isInfeasible = false;

		Stats->bestLowerBound = Cplex.getBestObjValue();
		if (!Stats->TSPNodes.empty()) {
			auto it = std::minmax_element(Stats->TSPTime.begin(), Stats->TSPTime.end());
			Stats->minTSPtime = *it.first;
			Stats->maxTSPTime = *it.second;
			Stats->avgTSPTime = std::accumulate(Stats->TSPTime.begin(), Stats->TSPTime.end(), 0.0) / Stats->TSPTime.size();
		}
		Stats->TSPSolved = Stats->TSPNodes.size();
		AddResultsToFile();
	} catch (IloException& e) {
	std::cerr << "Concert Exception: " << e << std::endl;
	}
	
}





int MFGVRP_Solver::NonDominatedChargingPaths() {
	//std::vector<std::vector<std::vector<int>>> R(n+1, std::vector<std::vector<int>>(n+1, std::vector<int>(1,NULL))); // Non-dominated charging paths.
	bool NonDominated;
	double MinDist;
	int MinChar;
	int numY = 0;
	int Depot = 0;
	for (IloInt i = 0; i < n + 1; i++)
	{
		Rj.clear();
		for (IloInt j = 0; j < n + 1; j++)
		{
			Rij.clear();
			if (i == 0 || j == 0) Depot = 1; else Depot = 0;
			for (IloInt r = n + 1 + Depot; r < n + f + 1; r++)
			{
				NonDominated = true;

				for (IloInt rr = n + 1+Depot; rr < n + f + 1; rr++)
				{
					if ((c[r][i][0] >= c[rr][i][0] && c[r][j][0] >= c[rr][j][0] && r != rr)) { NonDominated = false; break; }
				}
				if (NonDominated) { Rij.push_back(r); numY++; }
			}

			Rj.push_back(Rij);
		}
		MinDist = 9999;
		for (IloInt r = n + 1; r < n + f + 1; r++)
		{
			if (MinDist > c[r][i][0]) {
				MinDist = c[r][i][0];
				MinChar = r;
			}
		}
		if (i == 0)MinCharDist.push_back(0);
		else MinCharDist.push_back(MinDist); 
		R.push_back(Rj);
	}


	return numY;
}


void MFGVRP_Solver::AddResultsToFile() {
	std::ofstream outfile("results.txt", std::ios_base::app);
	double ICEV_Costs=0;
	double EV_Cost = 0;
	std::string name = Stats->instanceName;
	name.erase(name.find("."), 4);
	std::size_t PathID = name.find("\\");
	if (PathID != std::string::npos) name.erase(0, PathID + 1);
	std::ifstream out("results.txt", std::ios::app);
	outfile << "Name" << "\t" << "n" << "\t" << "r" << "\t" << "EVs" <<
		"\t" << "ICEVs" << "\t" << "Upper bound" << "\t" << "Lower bound" << "\t" << "Time" <<
		"\t" << "Optimal" << "\t" << "EV distance" << "\t" << "ICEV distance" << "\t"<< "GZ type" << "\t" << "GZ penalty" << "\t" << "Infeasible" << "\t" << "LB root" <<
		"\t" << "Time root" << "\t" << "Cap root" << "\t" <<
		"FCI root" << "\t" << "MS root" << "\t" << "Comb root" <<
		"\t" << "NCP root" << "\t" << "NCS root" << "\t" << "CP root"
		<< "\t" << "FC root" << "\t" << "FCS root"
		<< "\t" << "FE root" << "\t" << "FSE root"
		<< "\t" << "ICEV_T" << "\t" << "EV_T" << "\t" << "time_set" << "\t" << "time_EV_set"
		<< "\t" << "Tree size" << "\t" << "Max depth"	<< "\t" << "UB depth" << "\t" << "Feas depth"
		<< "\t" << "CAP" << "\t" <<
		"FCI" << "\t" << "MS" << "\t" << "Comb" <<
		"\t" << "NCP" << "\t" << "NCS" << "\t" << "CP"
		<< "\t" << "FC" << "\t" << "FCS"
		<< "\t" << "FE" << "\t" << "FSE"  << "\t" << "ICEV_T" <<"\t" << "EV_T" <<"\t" << "time_set" << "\t" << "time_EV_set" << "\t" << "Num TSP"
		<< "\t" << "Avg TSP" << "\t" << "Max TSP" << "\t" << "Min TSP" "\t" << "Total TSP time"<< "\t" << "Total Cap time" << "\t" <<"Total enumeration time" <<"\t" << "Total Sep time" << "\t" << "Fetch data time" << "\t" << "Add cuts time" << "\t" << "Convert vector time" << "\t" << "Seperation Vector time"
		<< std::endl;




	if (Stats->isOptimal)
	{
		IloNum3DMatrix xTemp(env, n+1);
		IloNum3DMatrix yTemp(env, n+1);


		IloNumArray xTemp1d(env, xDummy.getSize());
		Cplex.getValues(xTemp1d, xDummy);

		IloNumArray yTemp1d(env, yDummy.getSize());
		Cplex.getValues(yTemp1d, yDummy);


		int cnt = 0;

		for (int i = 0; i < n+1; i++)
		{
			xTemp[i] = IloNum2DMatrix(env, n+1);
			for (int j = 0; j < n+1; j++)
			{
				xTemp[i][j] = IloNumArray(env, 2);
				for (int k = 0; k < 2; k++)
				{
					xTemp[i][j][k] = xTemp1d[cnt];
					cnt++;
				}
				ICEV_Costs += c[i][j][1] * xTemp[i][j][1];
				EV_Cost += c[i][j][0] * xTemp[i][j][0];
			}
		}

		xTemp1d.end();
		cnt = 0;

		for (int i = 0; i < n+1; i++)
		{
			yTemp[i] = IloNum2DMatrix(env, n+1);
			for (int j = 0; j < n+1; j++)
			{
				yTemp[i][j] = IloNumArray(env, R[i][j].size());
				for (int r = 0; r < R[i][j].size(); r++)
				{
					yTemp[i][j][r] = yTemp1d[cnt];
					EV_Cost += (c[i][R[i][j][r]][0] + c[R[i][j][r]][j][0]) * yTemp[i][j][r];
					cnt++;
				}
			}
		}
		yTemp1d.end();


		std::ofstream Solfile("solutions/" + name + "_sol.json");
		int CurNode = 0;
		int NoOfRoutes = 0;
		Solfile << "{" << std::endl;
		Solfile << "\t\"Instance\":\t\"" << name << "\"," << std::endl;
		Solfile << "\t\"n\":\t" << n << "," << std::endl;
		Solfile << "\t\"r\":\t" << f << "," << std::endl;
		Solfile << "\t\"ObjVal\":\t" << Stats->bestUpperBound << "," << std::endl;
		Solfile << "\t\"EV_Cost\":\t" << EV_Cost << "," << std::endl;
		Solfile << "\t\"ICEV_Cost\":\t" << ICEV_Costs << "," << std::endl;
		Solfile << "\t\"noOfRoutesEVRoutes\":\t" << Stats->EVs << "," << std::endl;
		Solfile << "\t\"noOfRoutesICEVRoutes\":\t" << Stats->ICEVs << "," << std::endl;
		Solfile << "\t\"GZType\":\t" << GZType << "," << std::endl;
		Solfile << "\t\"GZPenalty\":\t" << GZPenalty << "," << std::endl;
		Solfile << "\t\"coordinates\":\t[";
		for (int i = 0; i < n+f; i++)
		{
			Solfile << "[" << xCoord[i] << "," << yCoord[i] << "],";
		}
		Solfile << "[" << xCoord[n+f] << "," << yCoord[n+f] << "]]," << std::endl;;
		if (G > 0) {
			Solfile << "\t\"GreenZones\":\t[";
			for (int i = 0; i < GreenZones.size() - 1; i++)
			{
				Solfile << "[" << GreenZones[i][0] << "," << GreenZones[i][1] << "," << GreenZones[i][2] << "," << GreenZones[i][3] << "],";
			}
			Solfile << "[" << GreenZones[GreenZones.size() - 1][0] << "," << GreenZones[GreenZones.size() - 1][1] << "," << GreenZones[GreenZones.size() - 1][2] << "," << GreenZones[GreenZones.size() - 1][3] << "]]," << std::endl;;
		}
		Solfile << "\t\"EVRoutes\":\t[";
		for (int i = 0; i < n + 1; i++)
		{
			if (xTemp[CurNode][i][0] > 0.9) {
				Solfile << "[" << CurNode << ", " << i;
				CurNode = i;
				for (int j = 0; j < n + 1; j++)
				{
					if (xTemp[CurNode][j][0] > 0.9) {
						Solfile << ", " << j;
						CurNode = j;
						j = -1;
					}
					if (j != -1)
					{
						for (int r = 0; r < R[CurNode][j].size(); r++)
						{
							if (yTemp[CurNode][j][r] > 0.9)
							{
								Solfile << ", " << R[CurNode][j][r] << ", " << j;
								CurNode = j;
								j = -1;
								break;
							}
						}
					}

					if (CurNode == 0)
					{
						NoOfRoutes++;
						if (NoOfRoutes < std::round(Stats->EVs)) {
							Solfile << "], ";
							break;
						}
						else
						{
							Solfile << "]";
							break;
						}


					}
				}

			}
			for (int r1 = 0; r1 < R[CurNode][i].size(); r1++)
			{
				if (yTemp[CurNode][i][r1] > 0.9) {
					Solfile << "[" << CurNode << ", " << R[CurNode][i][r1] << ", " << i;
					CurNode = i;
					for (int j = 0; j < n + 1; j++)
					{
						if (xTemp[CurNode][j][0] > 0.9) {
							Solfile << ", " << j;
							CurNode = j;
							j = -1;
						}
						if (j != -1) {
							for (int r = 0; r < R[CurNode][j].size(); r++)
							{
								if (yTemp[CurNode][j][r] > 0.9)
								{
									Solfile << ", " << R[CurNode][j][r] << ", " << j;
									CurNode = j;
									j = -1;
									break;
								}
							}
						}

						if (CurNode == 0)
						{
							NoOfRoutes++;
							if (NoOfRoutes < std::round(Stats->EVs)) {
								Solfile << "], ";
								break;
							}
							else
							{
								Solfile << "]";
								break;
							}
						}
					}
				}
			}



		}
		Solfile << "]," << std::endl;
		Solfile << "\t\"ICEVRoutes\":\t[";
		NoOfRoutes = 0;
		for (int i = 0; i < n + 1; i++)
		{
			if (xTemp[CurNode][i][1] > 0.9) {
				Solfile << "[" << CurNode << ", " << i;
				CurNode = i;
				for (int j = 0; j < n + 1; j++)
				{
					if (xTemp[CurNode][j][1] > 0.9) {
						Solfile << ", " << j;
						CurNode = j;
						j = -1;
					}

					if (CurNode == 0)
					{
						NoOfRoutes++;
						if (NoOfRoutes < std::round(Stats->ICEVs)) {
							Solfile << "], ";
							break;
						}
						else
						{
							Solfile << "]";
							break;
						}


					}
				}

			}
		}
		Solfile << "]" << std::endl;
		Solfile << "}";
		Solfile.close();

		for (int i = 0; i < n+1; i++)
		{
			for (int j = 0; j < n+1; j++)
			{
				xTemp[i][j].end();
			}
			xTemp[i].end();
		}
		xTemp.end();

		for (int i = 0; i < n+1; i++)
		{
			for (int j = 0; j < n+1; j++)
			{
				yTemp[i][j].end();
			}
			yTemp[i].end();
		}
		yTemp.end();
	}

	outfile << name << "\t" << Stats->n << "\t" << Stats->r << "\t" << Stats->EVs <<
		"\t" << Stats->ICEVs << "\t" << Stats->bestUpperBound << "\t" << Stats->bestLowerBound << "\t" << Stats->time <<
		"\t" << Stats->isOptimal << "\t" << EV_Cost << "\t" << ICEV_Costs << "\t" << GZType << "\t" << GZPenalty << "\t" << Stats->isInfeasible << "\t" << Stats->LBAtRootNode <<
		"\t" << Stats->timeRootNode << "\t" << Stats->numberOfCutsAtRootNode->CapCuts << "\t" <<
		Stats->numberOfCutsAtRootNode->FCI << "\t" << Stats->numberOfCutsAtRootNode->MultiStar << "\t" << Stats->numberOfCutsAtRootNode->Comb <<
		"\t" << Stats->numberOfCutsAtRootNode->NoChargePath << "\t" << Stats->numberOfCutsAtRootNode->NoChargeSet << "\t" << Stats->numberOfCutsAtRootNode->ChargerPath
		<< "\t" << Stats->numberOfCutsAtRootNode->FixedChargers << "\t" << Stats->numberOfCutsAtRootNode->FixedChargerSingle
		<< "\t" << Stats->numberOfCutsAtRootNode->FixedEdges << "\t" << Stats->numberOfCutsAtRootNode->FixedSingleEdge
		<< "\t" << Stats->numberOfCutsAtRootNode->TimeInfeasibleICEVPath << "\t" << Stats->numberOfCutsAtRootNode->TimeInfeasibleEVPath << "\t" << Stats->numberOfCutsAtRootNode->TimeInfeasibleSet << "\t" << Stats->numberOfCutsAtRootNode->TimeInfeasibleEVSet
		<< "\t" << Stats->treeSize << "\t" << Stats->treeDepth << "\t" << Stats->treeUBDepth << "\t" << Stats->treeFeasDepth
		<< "\t" << Stats->totalNumberOfCuts->CapCuts << "\t" <<
		Stats->totalNumberOfCuts->FCI << "\t" << Stats->totalNumberOfCuts->MultiStar << "\t" << Stats->totalNumberOfCuts->Comb <<
		"\t" << Stats->totalNumberOfCuts->NoChargePath << "\t" << Stats->totalNumberOfCuts->NoChargeSet << "\t" << Stats->totalNumberOfCuts->ChargerPath
		<< "\t" << Stats->totalNumberOfCuts->FixedChargers << "\t" << Stats->totalNumberOfCuts->FixedChargerSingle
		<< "\t" << Stats->totalNumberOfCuts->FixedEdges << "\t" << Stats->totalNumberOfCuts->FixedSingleEdge << "\t" << Stats->totalNumberOfCuts->TimeInfeasibleICEVPath << "\t" << Stats->totalNumberOfCuts->TimeInfeasibleEVPath << "\t" << Stats->totalNumberOfCuts->TimeInfeasibleSet << "\t" << Stats->totalNumberOfCuts->TimeInfeasibleEVSet << "\t" << Stats->TSPSolved
		<< "\t" << Stats->avgTSPTime << "\t" << Stats->maxTSPTime << "\t" << Stats->minTSPtime << "\t" << Stats->TSPSolved * Stats->avgTSPTime << "\t" << Stats->CapSepTime << "\t" << Stats->EnumerationTime - Stats->TSPSolved * Stats->avgTSPTime << "\t" << Stats->TotalSepTime << "\t" << Stats->TotalLoadSepData << "\t" << Stats->AddCutsTime << "\t" << Stats->ConvertVals << "\t" << Stats->SetupSepVectors
		<< std::endl;

	outfile.close();

}

void MFGVRP_Solver::SeparateCuts(UsedEdges* edges, CnstrMgrPointer MyCuts) {

	SeparateEnergyCuts(&(edges->DirectConnections), &(edges->ChargerConnections), MyCuts);

	if (MyCuts->Size == 0)
		SeparateDurationCuts(&(edges->ICEVConnections), &(edges->DirectConnections), &(edges->ChargerConnections), &(edges->Consolidated), MyCuts);

}

void MFGVRP_Solver::SeparateDurationCuts(std::vector<connections>* ICEV, std::vector<connections>* Direct, std::vector<connections>* Charger, std::vector<connections>* Consolidated, CnstrMgrPointer MyCuts) {

	double setSlack; double setSlackEV;
	double tournamentSlack; double tournamentSlackEV;
	double energyConsumption;
	double MinCharDuration;
	double TravelTime;
	double ServiceTime;
	int nChargersIncluded = 0;
	bool ChargerIncluded;
	std::vector<double> sumOfOutgoingArcs; std::vector<double> sumOfOutgoingArcsEV;
	std::vector<double> sumOfIngoingArcs; std::vector<double> sumOfIngoingArcsEV;


	std::vector<int>path;
	std::vector<std::vector<int>> Indices(3, std::vector<int>(n + 1, 0));
	double TSPval = 0;
	int curNode;
	int newNode = 0;
	int newCharger;
	int VehicleID;

	//Start with EVs:
	for (int i = 1; i <= n; i++)
	{

		path.clear();
		path.push_back(i);
		nChargersIncluded = 0;
		ChargerIncluded = false;
		sumOfIngoingArcs.clear();
		sumOfOutgoingArcs.clear();
		curNode = i;
		newCharger = 0;
		setSlack = 0; tournamentSlack = 0; setSlackEV = 0; tournamentSlackEV = 0;
		energyConsumption = 0; ServiceTime = s[i]; TravelTime = 0;
		for (int v = 0; v < 3; v++)
			Indices[v][i] = 0;

		do
		{
			while (Consolidated->at(curNode).succ.size() > Indices[0][curNode] || Direct->at(curNode).succ.size() > Indices[1][curNode] || Charger->at(curNode).succ.size() > Indices[2][curNode])
			{
				VehicleID = -1;
				if (Consolidated->at(curNode).succ.size() > Indices[0][curNode])
					if (std::find(path.begin(), path.end(), Consolidated->at(curNode).succ[Indices[0][curNode]]) == path.end() && Consolidated->at(curNode).succ[Indices[0][curNode]] > 0)
					{
						newNode = Consolidated->at(curNode).succ[Indices[0][curNode]];
						VehicleID = 0;
					}

				if (VehicleID == -1 && Direct->at(curNode).succ.size() > Indices[1][curNode])
					if (std::find(path.begin(), path.end(), Direct->at(curNode).succ[Indices[1][curNode]]) == path.end() && Direct->at(curNode).succ[Indices[1][curNode]] > 0)
					{
						newNode = Direct->at(curNode).succ[Indices[1][curNode]];
						VehicleID = 1;
					}

				if (VehicleID == -1 && Charger->at(curNode).succ.size() > Indices[2][curNode])
					if (std::find(path.begin(), path.end(), Charger->at(curNode).succ[Indices[2][curNode]]) == path.end() && Charger->at(curNode).succ[Indices[2][curNode]] > 0)
					{
						newNode = Charger->at(curNode).succ[Indices[2][curNode]];
						newCharger = Charger->at(curNode).succCharger[Indices[2][curNode]];
						VehicleID = 2;
						
					}
				if (VehicleID != -1)
				{
					sumOfIngoingArcs.push_back(0); sumOfOutgoingArcs.push_back(0); sumOfOutgoingArcsEV.push_back(0); sumOfIngoingArcsEV.push_back(0);
					//Find the sum of ingoing and outgoing arcs between path and new node:
					for (int v = 0; v < Consolidated->at(newNode).pre.size(); v++)
						if (std::find(path.begin(), path.end(), Consolidated->at(newNode).pre[v]) != path.end())
							sumOfIngoingArcs.back() += Consolidated->at(newNode).preVals[v];

					for (int v = 0; v < Direct->at(newNode).pre.size(); v++)
						if (std::find(path.begin(), path.end(), Direct->at(newNode).pre[v]) != path.end())
							sumOfIngoingArcsEV.back() += Direct->at(newNode).preVals[v];

					for (int v = 0; v < Charger->at(newNode).pre.size(); v++)
						if (std::find(path.begin(), path.end(), Charger->at(newNode).pre[v]) != path.end())
							sumOfIngoingArcsEV.back() += Charger->at(newNode).preVals[v];

					if ((tournamentSlack + 1 - sumOfIngoingArcs.back() < 0.999 && nChargersIncluded==0 && VehicleID!=2) || tournamentSlackEV + 1 - sumOfIngoingArcsEV.back()<0.999) {
						if (VehicleID <= 1)
						{
							path.push_back(newNode);
							if (nChargersIncluded==0 && c[0][i][0] + TravelTime + ServiceTime + c[curNode][newNode][0] + s[newNode] + c[newNode][0][0] > T) //&& tournamentSlack + 1 - sumOfIngoingArcs.back()<0.999)
							{
								// Add tournament constraint for ICEV:
								AddDurationTour(&path, 2, MyCuts);
								sumOfIngoingArcs.pop_back(); sumOfOutgoingArcs.pop_back(); sumOfIngoingArcsEV.pop_back(); sumOfOutgoingArcsEV.pop_back();
								path.pop_back();
								goto next_node;
							}
							else if (c[0][i][0] + TravelTime + ServiceTime + c[curNode][newNode][0] + s[newNode] + c[newNode][0][0] + std::max(0.0, (energyConsumption + c[curNode][newNode][0] + c[newNode][0][0]) * r) > T && tournamentSlackEV + 1 - sumOfIngoingArcsEV.back() < 0.999)
							{
								// Add tournament constraint for EV:
								AddDurationTour(&path, 0, MyCuts);
								sumOfIngoingArcs.pop_back(); sumOfOutgoingArcs.pop_back(); sumOfIngoingArcsEV.pop_back(); sumOfOutgoingArcsEV.pop_back();
								path.pop_back();
								goto next_node;
							}
							path.pop_back();
							// extend path
						}
						else if (tournamentSlackEV + 1 - sumOfIngoingArcsEV.back() - sumOfOutgoingArcsEV.back() < 0.999)
						{
							if (c[0][i][0] + TravelTime + ServiceTime + c[curNode][newCharger][0] + c[newCharger][newNode][0] + s[newNode] + c[newNode][0][0] + std::max(0.0, (energyConsumption + c[curNode][newCharger][0] + c[newCharger][newNode][0] + c[newNode][0][0] - B) * r) > T)
							{
								path.push_back(newCharger); path.push_back(newNode);
								AddDurationTour(&path, 0, MyCuts);
								sumOfIngoingArcs.pop_back(); sumOfOutgoingArcs.pop_back(); sumOfIngoingArcsEV.pop_back(); sumOfOutgoingArcsEV.pop_back();
								path.pop_back(); path.pop_back();
								goto next_node;

							}

						}
						tournamentSlack += 1 - sumOfIngoingArcs.back();
						tournamentSlackEV += 1 - sumOfIngoingArcsEV.back();

						if (VehicleID == 2)
						{
							path.push_back(newCharger); path.push_back(newNode);
							energyConsumption += c[curNode][newCharger][0] + c[newCharger][newNode][0];
							TravelTime += c[curNode][newCharger][0] + c[newCharger][newNode][0];
							nChargersIncluded++;
						}
						else
						{
							path.push_back(newNode);
							energyConsumption += c[curNode][newNode][0];
							TravelTime += c[curNode][newNode][0];

						}
						ServiceTime += s[newNode];

						Indices[VehicleID][newNode] = -1;

						Indices[VehicleID][curNode]++;


						curNode = newNode;

					}
					else {
						sumOfIngoingArcs.pop_back(); sumOfOutgoingArcs.pop_back(); sumOfIngoingArcsEV.pop_back(); sumOfOutgoingArcsEV.pop_back();
					}
				}

			next_node:
				if (VehicleID == -1)
					for (int k = 0; k < 3; k++)
						Indices[k][curNode]++;
				else
					Indices[VehicleID][curNode]++;
			}
			if (path.size() > 1) {
				if (path[path.size() - 2] > n || path[path.size() - 2] == 0) { //Handle back tracking differently if the last addition included a charger


					energyConsumption -= c[path[path.size() - 3]][path[path.size() - 2]][0] + c[path[path.size() - 2]][path.back()][0];
					TravelTime -= c[path[path.size() - 3]][path[path.size() - 2]][0] + c[path[path.size() - 2]][path.back()][0];
					ServiceTime -= s[path.back()];
					path.pop_back();
					nChargersIncluded--;

				}
				else
				{
					energyConsumption -= c[path[path.size() - 2]][path.back()][0];
					TravelTime -= c[path[path.size() - 2]][path.back()][0];
					ServiceTime -= s[path.back()];
				}

				tournamentSlack += sumOfIngoingArcs.back() - 1;
				tournamentSlackEV += sumOfIngoingArcsEV.back() - 1;
				sumOfIngoingArcs.pop_back();
				sumOfOutgoingArcs.pop_back();
				sumOfIngoingArcsEV.pop_back();
				sumOfOutgoingArcsEV.pop_back();
				curNode = path[path.size() - 2];



			}
			else {
				energyConsumption = 0;
				TravelTime = 0;
				ServiceTime = 0;
				setSlack = 0;
				setSlackEV = 0;
				tournamentSlack = 0;
				tournamentSlackEV = 0;
				sumOfIngoingArcs.clear();
				sumOfOutgoingArcs.clear();
				sumOfIngoingArcsEV.clear();
				sumOfOutgoingArcsEV.clear();
			}

			path.pop_back();

		} while (path.size() >= 1);

	}




	//Then do ICEVs:
	for (int i = 1; i <= n; i++)
	{

		path.clear();
		path.push_back(i);
		sumOfIngoingArcs.clear();
		sumOfOutgoingArcs.clear();
		curNode = i;
		newCharger = 0;
		tournamentSlack = 0;
		ServiceTime = s[i]; TravelTime = 0;
		for (int v = 0; v < 3; v++)
			Indices[v][i] = 0;

		do
		{
			while (ICEV->at(curNode).succ.size() > Indices[0][curNode])
			{
				VehicleID = -1;
				if (ICEV->at(curNode).succ.size() > Indices[0][curNode])
					if (std::find(path.begin(), path.end(), ICEV->at(curNode).succ[Indices[0][curNode]]) == path.end() && ICEV->at(curNode).succ[Indices[0][curNode]] > 0)
					{
						newNode = ICEV->at(curNode).succ[Indices[0][curNode]];
						VehicleID = 0;
					}
				if (VehicleID != -1)
				{
					sumOfIngoingArcs.push_back(0); sumOfOutgoingArcs.push_back(0);
					//Find the sum of ingoing and outgoing arcs between path and new node:
					for (int v = 0; v < ICEV->at(newNode).pre.size(); v++)
						if (std::find(path.begin(), path.end(), ICEV->at(newNode).pre[v]) != path.end())
							sumOfIngoingArcs.back() += ICEV->at(newNode).preVals[v];


					if (tournamentSlack + 1 - sumOfIngoingArcs.back() < 0.999) {//- sumOfOutgoingArcs.back() < 0.999) {

						path.push_back(newNode);
						if (c[0][i][1] + TravelTime + ServiceTime + c[curNode][newNode][1] + s[newNode] + c[newNode][0][1] > T)
						{
							// Add tournament constraint for ICEV:
							AddDurationTour(&path, 1, MyCuts);
							sumOfIngoingArcs.pop_back(); sumOfOutgoingArcs.pop_back();
							path.pop_back();
							goto next_icev_node;
						}

						tournamentSlack += 1 - sumOfIngoingArcs.back();
						TravelTime += c[curNode][newNode][1];

						ServiceTime += s[newNode];

						Indices[VehicleID][newNode] = -1;

						Indices[VehicleID][curNode]++;


						curNode = newNode;

					}
					else {
						sumOfIngoingArcs.pop_back(); sumOfOutgoingArcs.pop_back();
					}
				}

			next_icev_node:
				if (VehicleID == -1)
					for (int k = 0; k < 3; k++)
						Indices[k][curNode]++;
				else
					Indices[VehicleID][curNode]++;
			}
			if (path.size() > 1) {
				TravelTime -= c[path[path.size() - 2]][path.back()][1];
				ServiceTime -= s[path.back()];
				tournamentSlack += sumOfIngoingArcs.back() - 1;
				sumOfIngoingArcs.pop_back();
				sumOfOutgoingArcs.pop_back();
				curNode = path[path.size() - 2];
			}
			else {
				energyConsumption = 0;
				TravelTime = 0;
				ServiceTime = 0;
				tournamentSlack = 0;
				sumOfIngoingArcs.clear();
				sumOfOutgoingArcs.clear();
			}

			path.pop_back();

		} while (path.size() >= 1);

	}



}

void MFGVRP_Solver::SeparateEnergyCuts(std::vector<connections>* Direct, std::vector<connections>* Charger, CnstrMgrPointer MyCuts) {

	std::vector<int>path;
	std::vector<std::pair<int, double>> IngoingChargersPath;
	std::vector<std::pair<int, double>> OutgoingChargersPath;
	std::vector<std::pair<int, double>> IngoingChargersSet;
	std::vector<std::pair<int, double>> OutgoingChargersSet;
	std::vector<std::pair<int, double>> IngoingDualChargers;
	std::vector<std::pair<int, double>> OutgoingDualChargers;
	double setSlack;
	double tournamentSlack;
	double dualPathSlack;
	double energyConsumption;
	std::vector<double> sumOfOutgoingArcs;
	std::vector<double> sumOfIngoingArcs;
	std::vector<double> sumOfDualArcs;
	int curNode;
	int newNode;
	bool CheckChargers = true;
	bool ChargerConnected;
	std::vector<int> DirectIndex(n + 1, 0);
	double TSPval = 0;
	//Start building paths/set from direct edges used in the current LP solution:
	for (int i = 1; i < n + 1; i++)
	{

		//Try each customer as the source node
		if (Direct->at(i).succ.size() == 0)
			continue;

		path.clear();
		path.push_back(i);
		sumOfIngoingArcs.clear();
		sumOfOutgoingArcs.clear();
		DirectIndex[i] = 0;
		curNode = i;
		setSlack = 0; tournamentSlack = 0; energyConsumption = 0; dualPathSlack = 0;
		do
		{
			while (Direct->at(curNode).succ.size() > DirectIndex[curNode])
			{
				//If the node we try to add is not the depot and not already present in the path:
				if (Direct->at(curNode).succ[DirectIndex[curNode]] != 0 && std::find(path.begin(), path.end(), Direct->at(curNode).succ[DirectIndex[curNode]]) == path.end())
				{

					newNode = Direct->at(curNode).succ[DirectIndex[curNode]];
					//DirectIndex[newNode] = 0;
					sumOfIngoingArcs.push_back(0); sumOfOutgoingArcs.push_back(0); sumOfDualArcs.push_back(0);
					//Find the sum of ingoing and outgoing arcs between path and new node:
					for (int v = 0; v < Direct->at(newNode).pre.size(); v++)
						if (std::find(path.begin(), path.end(), Direct->at(newNode).pre[v]) != path.end()) {
							sumOfIngoingArcs.back() += Direct->at(newNode).preVals[v];
							if (Direct->at(newNode).pre[v] == path.back())
								sumOfDualArcs.back() += Direct->at(newNode).preVals[v];
						}


					for (int v = 0; v < Direct->at(newNode).succ.size(); v++)
						if (std::find(path.begin(), path.end(), Direct->at(newNode).succ[v]) != path.end()) {
							sumOfOutgoingArcs.back() += Direct->at(newNode).succVals[v];
							if (Direct->at(newNode).succ[v] == path.back())
								sumOfDualArcs.back() += Direct->at(newNode).succVals[v];
						}

					//If the set slack of union of the nodes in the path and the new node allows us to find a violated inequality, we continue with the path
					if (setSlack + 1 - sumOfIngoingArcs.back() - sumOfOutgoingArcs.back() < 0.999)
					{
						if (energyConsumption + c[curNode][newNode][0] + MinCharDist[path[0]] + MinCharDist[newNode] > B)
						{
							path.push_back(newNode);

							if (dualPathSlack + 1 - sumOfDualArcs.back() < 0.999)
							{
								AddNoChargePath(&path, MyCuts);
								sumOfIngoingArcs.pop_back(); sumOfOutgoingArcs.pop_back(); sumOfDualArcs.pop_back();

								path.pop_back();
								goto next_node;
							}


							//If we can cut it off with a tournament constraint, we do that:
							if (tournamentSlack + 1 - sumOfIngoingArcs.back() < 0.999)
							{
								// ----- Add tournament constraint to LP -----
								AddNoChargeTour(&path, MyCuts);
								sumOfIngoingArcs.pop_back(); sumOfOutgoingArcs.pop_back(); sumOfDualArcs.pop_back();
								path.pop_back();
								goto next_node;
							}
							else //Check if we have a valid no charge set inequality:
							{
								// ----- Solve TSP to validate no charge set inequality

								TSPval = TSPConcorde(&path);
								if (TSPval > B)
								{
									AddNoChargeSet(&path, MyCuts);
									sumOfIngoingArcs.pop_back(); sumOfOutgoingArcs.pop_back(); sumOfDualArcs.pop_back();
									path.pop_back();
									goto next_node;
								}
							}
							path.pop_back();
						}
						// Check if we have a valid set/path when we add chargers


						if (CheckChargers)
						{
							//#### CHECK PATH FIRST TO CHARGERS BEFORE STARTING THE VALIDATION OF OTHER INEQUALITIES (TOURNEMENT)



							//Check if we have a violation

							// Get connections between chargers and path:
							path.push_back(newNode);
							IngoingChargersPath.clear(); OutgoingChargersPath.clear(); IngoingChargersSet.clear(); OutgoingChargersSet.clear();
							CheckChargerConnections(Charger, &(Direct->at(0)), &path, IngoingChargersPath, OutgoingChargersPath, IngoingChargersSet, OutgoingChargersSet);


							//Check ingoing tournament:
							for (int r1 = 0; r1 < IngoingChargersPath.size(); r1++)
							{
								if (tournamentSlack + 2 - sumOfIngoingArcs.back() - IngoingChargersPath[r1].second < 0.999 && energyConsumption + c[curNode][path.back()][0] + MinCharDist[path.back()] + c[IngoingChargersPath[r1].first][path[0]][0]>B)
								{
									// Add violated tournament constriant for r1
									AddChargerTour(&path, IngoingChargersPath[r1].first, -1, MyCuts);
									sumOfIngoingArcs.pop_back(); sumOfOutgoingArcs.pop_back(); sumOfDualArcs.pop_back();
									path.pop_back();
									goto next_node;
								}
							}

							//Check outgoing tournament:
							for (int r2 = 0; r2 < OutgoingChargersPath.size(); r2++)
							{
								if (tournamentSlack + 2 - sumOfOutgoingArcs.back() - OutgoingChargersPath[r2].second < 0.999 && energyConsumption + c[curNode][path.back()][0] + MinCharDist[path[0]] + c[path.back()][OutgoingChargersPath[r2].first][0]>B)
								{
									// Add violated tournament constriant for r2
									AddChargerTour(&path, -1, OutgoingChargersPath[r2].first, MyCuts);
									sumOfIngoingArcs.pop_back(); sumOfOutgoingArcs.pop_back(); sumOfDualArcs.pop_back();
									path.pop_back();
									goto next_node;
								}
							}

							//Check ingoing/outgoing pair for tournament:
							for (int r1 = 0; r1 < IngoingChargersPath.size(); r1++)
							{
								ChargerConnected = false;
								if (IngoingChargersPath[r1].first == 0) {
									if (std::find(Direct->at(0).succ.begin(), Direct->at(0).succ.end(), path[0]) != Direct->at(0).succ.end())
										ChargerConnected = true;
								}
								else if (std::find(Charger->at(path[0]).preCharger.begin(), Charger->at(path[0]).preCharger.end(), IngoingChargersPath[r1].first) != Charger->at(path[0]).preCharger.end())
									ChargerConnected = true;

								if (!ChargerConnected) continue;

								for (int r2 = 0; r2 < OutgoingChargersPath.size(); r2++)
								{
									ChargerConnected = false;
									if (OutgoingChargersPath[r2].first == 0) {
										if (std::find(Direct->at(0).pre.begin(), Direct->at(0).pre.end(), path.back()) != Direct->at(0).pre.end())
											ChargerConnected = true;
									}
									else if (std::find(Charger->at(path.back()).succCharger.begin(), Charger->at(path.back()).succCharger.end(), OutgoingChargersPath[r2].first) != Charger->at(path.back()).succCharger.end())
										ChargerConnected = true;

									if (!ChargerConnected) continue;

									if (tournamentSlack + 3 - sumOfOutgoingArcs.back() - IngoingChargersPath[r1].second - OutgoingChargersPath[r2].second < 0.999 && energyConsumption + c[curNode][path.back()][0] + c[IngoingChargersPath[r1].first][path[0]][0] + c[path.back()][OutgoingChargersSet[r2].first][0]>B)
									{
										// Add violated tournament constriant for r1 r2 pair
										AddChargerTour(&path, IngoingChargersPath[r1].first, OutgoingChargersPath[r2].first, MyCuts);
										sumOfIngoingArcs.pop_back(); sumOfOutgoingArcs.pop_back(); sumOfDualArcs.pop_back();
										path.pop_back();
										goto next_node;
									}
								}
							}


							//Check single charger edges and single chargers:
							//See if charger violates the path
							for (int r1 = 0; r1 < IngoingChargersSet.size(); r1++)
							{
								if ((path.size()) * 2 - (setSlack + 1 - sumOfIngoingArcs.back()) * 2 - 2 + IngoingChargersSet[r1].second > path.size() * 2 - 1.999 && energyConsumption + c[curNode][path.back()][0] + MinCharDist[path.back()] + c[IngoingChargersSet[r1].first][path[0]][0] > B)
								{
									//We have a violated set. Now check if we can validate it
									TSPval = TSPConcorde(&path,0, false, false,false, IngoingChargersSet[r1].first, -1);
									if (TSPval > B) {
										AddFixedCharger(&path, IngoingChargersSet[r1].first, -1, MyCuts);
										sumOfIngoingArcs.pop_back(); sumOfOutgoingArcs.pop_back(); sumOfDualArcs.pop_back();
										path.pop_back();
										goto next_node;
									}
								}
							}
							for (int r2 = 0; r2 < OutgoingChargersSet.size(); r2++)
							{
								if ((path.size()) * 2 - (setSlack + 1 - sumOfOutgoingArcs.back()) * 2 - 2 + OutgoingChargersSet[r2].second > path.size() * 2 - 1.999 && energyConsumption + c[curNode][path.back()][0] + MinCharDist[path[0]] + c[path.back()][OutgoingChargersSet[r2].first][0] > B)
								{

									////We have a violated set. Now check if we can validate it
									TSPval = TSPConcorde(&path,0, false, false,false, -1, OutgoingChargersSet[r2].first);
									if (TSPval > B) {
										AddFixedCharger(&path, -1, OutgoingChargersSet[r2].first, MyCuts);
										sumOfIngoingArcs.pop_back(); sumOfOutgoingArcs.pop_back(); sumOfDualArcs.pop_back();
										path.pop_back();
										goto next_node;
									}
								}
							}
							//Try to fixate edge
							for (int r1 = 0; r1 < IngoingChargersPath.size(); r1++)
							{
								if (path.size() - 1 - (setSlack + 1 - sumOfIngoingArcs.back()) > path.size() - 0.999 - IngoingChargersPath[r1].second && energyConsumption + c[curNode][path.back()][0] + MinCharDist[path.back()] + c[IngoingChargersPath[r1].first][path[0]][0] > B)
								{
									////We have a violated set. Now check if we can validate it
									TSPval = TSPConcorde(&path,0, false, false,true, IngoingChargersPath[r1].first, -1);
									if (TSPval > B) {
										AddFixedChargingEdge(&path, IngoingChargersPath[r1].first, -1, MyCuts);
										sumOfIngoingArcs.pop_back(); sumOfOutgoingArcs.pop_back(); sumOfDualArcs.pop_back();
										path.pop_back();
										goto next_node;
									}
								}
							}

							for (int r2 = 0; r2 < OutgoingChargersPath.size(); r2++)
							{
								if (path.size() - 1 - (setSlack + 1 - sumOfOutgoingArcs.back()) > path.size() - 0.999 - OutgoingChargersPath[r2].second && energyConsumption + c[curNode][path.back()][0] + MinCharDist[path[0]] + c[path.back()][OutgoingChargersPath[r2].first][0] > B)
								{
									//We have a violated set. Now check if we can validate it
									TSPval = TSPConcorde(&path, 0, false, false,true, -1,OutgoingChargersPath[r2].first);
									if (TSPval > B) {
										AddFixedChargingEdge(&path, -1, OutgoingChargersPath[r2].first, MyCuts);
										sumOfIngoingArcs.pop_back(); sumOfOutgoingArcs.pop_back(); sumOfDualArcs.pop_back();
										path.pop_back();
										goto next_node;
									}
								}
							}


							//Check pair of charger edges and pair of chargers:
							//See if a set of chargers violates the path
							for (int r1 = 0; r1 < IngoingChargersSet.size(); r1++)
							{
								//Make a condition here, if ingoing charger value is too small for finding a violation
								for (int r2 = 0; r2 < OutgoingChargersSet.size(); r2++)
								{
									if ((path.size()) * 3 - (setSlack + 1 - sumOfIngoingArcs.back() - sumOfOutgoingArcs.back()) * 3 - 3 + IngoingChargersSet[r1].second + OutgoingChargersSet[r2].second > path.size() * 3 - 1.999 && energyConsumption + c[curNode][path.back()][0] + c[IngoingChargersSet[r1].first][path[0]][0] + c[path.back()][OutgoingChargersSet[r2].first][0] > B)
									{
										////We have a violated set. Now check if we can validate it
										TSPval = TSPConcorde(&path,0, false, false,false, IngoingChargersSet[r1].first, OutgoingChargersSet[r2].first);
										if (TSPval > B) {
											AddFixedCharger(&path, IngoingChargersSet[r1].first, OutgoingChargersSet[r2].first, MyCuts);
											sumOfIngoingArcs.pop_back(); sumOfOutgoingArcs.pop_back(); sumOfDualArcs.pop_back();
											path.pop_back();
											goto next_node;
										}
									}
								}
							}


							//Try to fixate edges
							for (int r1 = 0; r1 < IngoingChargersPath.size(); r1++)
							{
								//Make a condition here, if ingoing charger value is too small for finding a violation
								for (int r2 = 0; r2 < OutgoingChargersPath.size(); r2++)
								{
									if (path.size() - 1 - (setSlack + 1 - sumOfIngoingArcs.back() - sumOfOutgoingArcs.back()) > path.size() - IngoingChargersPath[r1].second - OutgoingChargersPath[r2].second + 0.001 && energyConsumption + c[curNode][path.back()][0] + c[IngoingChargersPath[r1].first][path[0]][0] + c[path.back()][OutgoingChargersPath[r2].first][0] > B)
									{
										//We have a violated set. Now check if we can validate it
										TSPval = TSPConcorde(&path,0, false, false,true, IngoingChargersPath[r1].first, OutgoingChargersPath[r2].first);
										if (TSPval > B) {
											AddFixedChargingEdge(&path, IngoingChargersPath[r1].first, OutgoingChargersPath[r2].first, MyCuts);
											sumOfIngoingArcs.pop_back(); sumOfOutgoingArcs.pop_back(); sumOfDualArcs.pop_back();
											path.pop_back();
											goto next_node;
										}
									}
								}
							}

							path.pop_back();
						}


						// else just extend path:
						tournamentSlack += 1 - sumOfIngoingArcs.back();
						setSlack += 1 - sumOfIngoingArcs.back() - sumOfOutgoingArcs.back();
						dualPathSlack += 1 - sumOfDualArcs.back();
						energyConsumption += c[curNode][newNode][0];
						path.push_back(newNode);
						DirectIndex[newNode] = -1;
						DirectIndex[curNode]++;

						curNode = newNode;


					}
					else {
						sumOfIngoingArcs.pop_back(); sumOfOutgoingArcs.pop_back(); sumOfDualArcs.pop_back();
					}



				}
			next_node:
				DirectIndex[curNode]++;
			}
			//printf("\n\n");
			//for (int t = 0; t < path.size(); t++)
			//{
			//	printf("%d ", path[t]);
			//}

			if (path.size() > 1) {
				tournamentSlack += sumOfIngoingArcs.back() - 1;
				setSlack += sumOfIngoingArcs.back() + sumOfOutgoingArcs.back() - 1;
				dualPathSlack += sumOfDualArcs.back() - 1;
				energyConsumption -= c[path[path.size() - 2]][path.back()][0];
				curNode = path[path.size() - 2];
				sumOfIngoingArcs.pop_back();
				sumOfOutgoingArcs.pop_back();
				sumOfDualArcs.pop_back();

			}
			else {
				energyConsumption = 0;
				setSlack = 0;
				tournamentSlack = 0;
				sumOfIngoingArcs.clear();
				sumOfOutgoingArcs.clear();
				sumOfDualArcs.clear();
			}

			path.pop_back();


		} while (path.size() >= 1);
	}
	return;
}


void MFGVRP_Solver::CheckChargerConnections(std::vector<connections>* Charger, connections* Depot, std::vector<int>* path, std::vector<std::pair<int, double>>& EligibleIngoingChargers, std::vector<std::pair<int, double>>& EligibleOutgoingChargers, std::vector<std::pair<int, double>>& EligibleIngoingChargersSet, std::vector<std::pair<int, double>>& EligibleOutgoingChargersSet) {




	std::vector<double> IngoingChargers(f + 1, 0.0);
	std::vector<double> OutgoingChargers(f + 1, 0.0);
	std::vector<double> IngoingChargersSet(f + 1, 0.0);
	std::vector<double> OutgoingChargersSet(f + 1, 0.0);

	//Check connections to/from depot:

		//First specific edges:
	for (int v = 0; v < Depot->succ.size(); v++)
		if (Depot->succ[v] == path->at(0))
			IngoingChargers[0] += Depot->succVals[v];

	for (int v = 0; v < Depot->pre.size(); v++)
		if (Depot->pre[v] == path->back())
			OutgoingChargers[0] += Depot->preVals[v];


	//Then any connection between depot and the path:
	for (int u = 0; u < path->size(); u++)
	{
		for (int v = 0; v < Depot->succ.size(); v++)
			if (Depot->succ[v] == path->at(u))
				IngoingChargersSet[0] += Depot->succVals[v];
		for (int v = 0; v < Depot->pre.size(); v++)
			if (Depot->pre[v] == path->at(u))
				OutgoingChargersSet[0] += Depot->preVals[v];
	}

	//Check connections to/from chargers:
		/*Charger->at(path->at(0))*/
		//First specific edges
	if (Charger->at(path->at(0)).preCharger.size() > 0)
		for (int v = 0; v < Charger->at(path->at(0)).preCharger.size(); v++)
			IngoingChargers[Charger->at(path->at(0)).preCharger[v] - n] += Charger->at(path->at(0)).preVals[v];

	if (Charger->at(path->back()).succCharger.size() > 0)
		for (int v = 0; v < Charger->at(path->back()).succCharger.size(); v++)
			OutgoingChargers[Charger->at(path->back()).succCharger[v] - n] += Charger->at(path->back()).succVals[v];



	//Then check any connections between depot and path:
	for (int u = 0; u < path->size(); u++)
	{
		if (Charger->at(path->at(u)).preCharger.size() > 0)
			for (int v = 0; v < Charger->at(path->at(u)).preCharger.size(); v++)
				IngoingChargersSet[Charger->at(path->at(u)).preCharger[v] - n] += Charger->at(path->at(u)).preVals[v];
		if (Charger->at(path->at(u)).succCharger.size() > 0)
			for (int v = 0; v < Charger->at(path->at(u)).succCharger.size(); v++)
				OutgoingChargersSet[Charger->at(path->at(u)).succCharger[v] - n] += Charger->at(path->at(u)).succVals[v];
	}


	//prepare result vectors


	for (int id = 0; id < f + 1; id++)
	{
		if (id == 0)
		{
			if (IngoingChargers[id] > 0.0)
				EligibleIngoingChargers.push_back(std::make_pair(id, IngoingChargers[id]));

			if (OutgoingChargers[id] > 0.0)
				EligibleOutgoingChargers.push_back(std::make_pair(id, OutgoingChargers[id]));

			if (IngoingChargersSet[id] > 0.0)
				EligibleIngoingChargersSet.push_back(std::make_pair(id, IngoingChargersSet[id]));

			if (OutgoingChargersSet[id] > 0.0)
				EligibleOutgoingChargersSet.push_back(std::make_pair(id, OutgoingChargersSet[id]));
		}
		else
		{
			if (IngoingChargers[id] > 0.0)
				EligibleIngoingChargers.push_back(std::make_pair(id + n, IngoingChargers[id]));

			if (OutgoingChargers[id] > 0.0)
				EligibleOutgoingChargers.push_back(std::make_pair(id + n, OutgoingChargers[id]));

			if (IngoingChargersSet[id] > 0.0)
				EligibleIngoingChargersSet.push_back(std::make_pair(id + n, IngoingChargersSet[id]));

			if (OutgoingChargersSet[id] > 0.0)
				EligibleOutgoingChargersSet.push_back(std::make_pair(id + n, OutgoingChargersSet[id]));
		}

	}
	return;
}

double MFGVRP_Solver::TSPConcorde(std::vector<int>* path,int VehicleIndex, bool Duration, bool symmetric, bool fixed, int charg_in, int charg_out) {


	std::vector<std::vector<double>> dist(path->size() + 1, std::vector<double>(path->size() + 1));
	if (Duration)
		setupTSPDurationMatrix(path, dist,VehicleIndex);
	else
		setupTSPEnergyMatrix(path, dist, symmetric, fixed, charg_in, charg_out);

	if (!symmetric)
	{
		double objVal;
		int nCount = (path->size() + 1) * 3;
		int M = 999;
		int ecount = (nCount * (nCount - 1)) / 2; //Number of edges
		int* elist = new int[ecount * 2]; //Array giving the ends of the edges (in pairs)
		int* elen = new int[ecount]; //Array giving the weights of the edges
		int edge = 0;
		int edgeWeight = 0;
		int nodetype_i = 0;
		int nodetype_j = 0;
		int org_i = 0;
		int org_j = 0;
		for (int i = 0; i < nCount; i++)
		{
			//printf("\n%d", i);
			for (int j = i + 1; j < nCount; j++)
			{
				nodetype_i = i % 3;
				nodetype_j = j % 3;
				org_i = i / 3;
				org_j = j / 3;
				if (nodetype_i == nodetype_j)
				{
					elist[edge] = i;
					elist[edge + 1] = j;
					elen[edgeWeight] = M * 1000;
					edgeWeight++;
					edge = edge + 2;
				}
				else if (i + 1 == j && org_i == org_j)
				{
					elist[edge] = i;
					elist[edge + 1] = j;
					elen[edgeWeight] = 0;
					edgeWeight++;
					edge = edge + 2;
				}
				else if (org_i != org_j)
				{
					if (nodetype_i == 0 && nodetype_j == 2)
					{
						elist[edge] = i;
						elist[edge + 1] = j;
						elen[edgeWeight] = dist[org_j][org_i] * 1000;
						edgeWeight++;
						edge = edge + 2;
					}
					else if (nodetype_i == 2 && nodetype_j == 0)
					{
						elist[edge] = i;
						elist[edge + 1] = j;
						elen[edgeWeight] = dist[org_i][org_j] * 1000;
						edgeWeight++;
						edge = edge + 2;
					}
					else
					{
						elist[edge] = i;
						elist[edge + 1] = j;
						elen[edgeWeight] = M * 1000;
						edgeWeight++;
						edge = edge + 2;
					}
				}
				else
				{
					elist[edge] = i;
					elist[edge + 1] = j;
					elen[edgeWeight] = M * 1000;
					edgeWeight++;
					edge = edge + 2;
				}
			}
		}

		//Run TSP
		//printf("\nstart condorde");
		solveTSPConcorde(elist, elen, nCount, ecount, objVal);
		//printf("\nend condorde");
		//Release memory
		//printf("\nStart releasning");
		delete[] elist;
		delete[] elen;
		//printf("\nEnd Releasing");
		return objVal / 1000;
	}
	return 0;

}

void MFGVRP_Solver::setupTSPEnergyMatrix(std::vector<int>* path, std::vector<std::vector<double>>& dist, bool symmetric, bool fixed, int charg_in, int charg_out) {

	int M = 999;
	for (int i = 0; i <= path->size(); i++)
	{
		for (int j = 0; j <= path->size(); j++)
		{
			if (i == j)
				continue;
			if (i == 0)
			{
				if (charg_in >= 0) {
					if (fixed)
					{
						if (j == 1)
						{
							dist[i][j] = c[charg_in][path->at(0)][0];
						}
						else
						{
							dist[i][j] = M;
						}
					}
					else
					{
						dist[i][j] = c[charg_in][path->at(j - 1)][0];
					}
				}
				else
				{
					dist[i][j] = MinCharDist[path->at(j - 1)];
				}
			}
			else if (j == 0)
			{
				if (charg_out >= 0) {
					if (fixed)
					{
						if (i == path->size())
						{
							dist[i][j] = c[path->back()][charg_out][0];
						}
						else
						{
							dist[i][j] = M;
						}
					}
					else
					{
						dist[i][j] = c[path->at(i - 1)][charg_out][0];
					}
				}
				else
				{
					dist[i][j] = MinCharDist[path->at(i - 1)];
				}
			}
			else
			{
				dist[i][j] = c[path->at(i - 1)][path->at(j - 1)][0];
			}
		}
	}
}

void MFGVRP_Solver::setupTSPDurationMatrix(std::vector<int>* path, std::vector<std::vector<double>>& dist,int VehicleIndex) {
	int M = 999;
	for (int i = 0; i <= path->size(); i++)
	{
		for (int j = 0; j <= path->size(); j++)
		{
			if (i == j)
				continue;
			if (i == 0)
			{
				dist[i][j] = c[0][path->at(j - 1)][VehicleIndex];
			}
			else if (j == 0)
			{
				dist[i][j] = c[path->at(i - 1)][0][VehicleIndex];
			}
			else
			{
				dist[i][j] = c[path->at(i - 1)][path->at(j - 1)][VehicleIndex];
			}


		}
	}
}
void MFGVRP_Solver::solveTSPConcorde(int* elist, int* elen, int nCities, int nEdges, double& optVal) {

	int rval = 0; //Concorde functions return 1 if something fails
	double* in_val = (double*)NULL; //Can be used to specify an initial upperbound (it can be NULL)
	double* timebound = (double*)NULL;; //Run time limit
	int success; //1 if the run finished normally, and set to 0 if the search was terminated early (by hitting some predefined limit) 
	int foundtour; //1 if a tour has been found (if success is 0, then it may not be the optimal tour)   
	int hit_timebound = 0; //1 if timebound was reached
	int* in_tour = (int*)NULL; //Gives a starting tour in node node node format (it can be NULL)
	int* out_tour = (int*)NULL; //Optimal tour (it can be NULL, if it is not NULL then it should point to an array of length at least ncount).  
	char* name = (char*)NULL; //Specifes a char string that will be used to name various files that are written during the branch and bound search
	static int silent = 1; //Suppress most output if set to a nonzero value
	CCrandstate rstate;
	int seed = rand();
	CCutil_sprand(seed, &rstate); //Initialize the portable random number generator


	out_tour = CC_SAFE_MALLOC(nCities, int);
	name = CCtsp_problabel("_");
	CCdatagroup dat;

	//Initialize a CCdatagroup
	CCutil_init_datagroup(&dat);

	//Convert a matrix of edge lengths to a CCdatagroup
	rval = CCutil_graph2dat_matrix(nCities, nEdges, elist, elen, 1, &dat);


	//Solves the TSP over the graph specified in the datagroup
	rval = CCtsp_solve_dat(1, nCities, &dat, in_tour, out_tour, in_val, &optVal, &success, &foundtour, name, timebound, &hit_timebound, silent, &rstate);
	CC_FREE(out_tour, int);
	CC_FREE(name, char);

}
void MFGVRP_Solver::FinalFeasibilityCheck(std::vector<connections>* Direct, std::vector<connections>* Charger, std::vector<connections>* ICEV, CnstrMgrPointer MyCuts) {
	std::vector<int> path;
	std::vector<int> PathFromCharger;
	double energyConsumption;
	double energyConsumptionSegment;
	double duration;
	double minCharTime = 0;
	bool ChargerInPath = false;
	int curNode;
	int newNode;
	int newCharger;
	int previousCharger;
	//Check ICEV Routes:
	for (int i = 0; i < ICEV->at(0).succ.size(); i++)
	{
		curNode = ICEV->at(0).succ[i];
		path.push_back(curNode);
		duration = c[0][curNode][1] + s[curNode];

		while (curNode != 0)
		{
			newNode = ICEV->at(curNode).succ[0];
			if (newNode != 0)
				path.push_back(newNode);

			duration += c[curNode][newNode][1] + s[newNode];
			curNode = newNode;
			if (duration > T)
			{
				AddDurationTour(&path, 1, MyCuts);
				break;
			}
		}
		path.clear();
	}
	//Check EV routes that start with a direct edge:
	for (int i = 0; i < Direct->at(0).succ.size(); i++)
	{
		path.clear();
		//PathFromCharger.push_back(0);
		previousCharger = 0;
		curNode = Direct->at(0).succ[i];
		path.push_back(curNode);
		PathFromCharger.push_back(curNode);
		energyConsumption = c[0][curNode][0];
		energyConsumptionSegment = c[0][curNode][0];
		duration = c[0][curNode][0] + s[curNode];
		ChargerInPath = false;
		while (curNode != 0)
		{
			if (Direct->at(curNode).succ.size() > 0)
			{
				newNode = Direct->at(curNode).succ[0];
				if (newNode != 0) {
					path.push_back(newNode);
					PathFromCharger.push_back(newNode);
				}
				else
					newCharger = 0;



				duration += c[curNode][newNode][0] + s[newNode];
				energyConsumption += c[curNode][newNode][0];
				energyConsumptionSegment += c[curNode][newNode][0];
				minCharTime = std::max(0.0, (energyConsumption - B) * r);
				curNode = newNode;

				if (duration + minCharTime > T)
				{
					if (ChargerInPath)
						AddDurationTour(&path, 0, MyCuts);
					else
						AddDurationTour(&path, 1, MyCuts);
					break;

				}
				if (energyConsumptionSegment > B)
				{
					if (newNode == 0)
						AddChargerTour(&PathFromCharger, previousCharger, newCharger, MyCuts);
					else
						AddChargerTour(&PathFromCharger, previousCharger, -1, MyCuts);
					break;
				}
			}
			else if (Charger->at(curNode).succ.size() > 0)
			{
				newNode = Charger->at(curNode).succ[0];
				newCharger = Charger->at(curNode).succCharger[0];
				energyConsumption += c[curNode][newCharger][0] + c[newCharger][newNode][0];
				energyConsumptionSegment += c[curNode][newCharger][0];
				duration += c[curNode][newCharger][0] + c[newCharger][newNode][0];
				path.push_back(newCharger);

				if (newNode != 0)
					path.push_back(newNode);
				curNode = newNode;
				minCharTime = std::max(0.0, (energyConsumption - B) * r);
				ChargerInPath = true;
				if (duration + minCharTime > T)
				{
					AddDurationTour(&path, 0, MyCuts);
					break;
				}
				if (energyConsumptionSegment > B)
				{
					AddChargerTour(&PathFromCharger, previousCharger, newCharger, MyCuts);
					break;
				}
				PathFromCharger.clear();
				previousCharger = newCharger;
				PathFromCharger.push_back(newNode);
				energyConsumptionSegment = c[newCharger][newNode][0];
			}
			else
				break;

		}
		path.clear();
		PathFromCharger.clear();
	}


	//Check EV routes that start with a charging edge:
	for (int i = 0; i < Charger->at(0).succ.size(); i++)
	{
		path.clear();
		PathFromCharger.clear();
		curNode = Charger->at(0).succ[i];
		newCharger = Charger->at(0).succCharger[i];
		path.push_back(newCharger);
		if (c[0][newCharger][0] > B)
		{
			AddChargerTour(&path, 0, newCharger, MyCuts);
			continue;
		}
		path.push_back(curNode);
		previousCharger = newCharger;
		PathFromCharger.push_back(curNode);
		energyConsumption = c[0][newCharger][0] + c[newCharger][curNode][0];
		energyConsumptionSegment = c[newCharger][curNode][0];
		duration = c[0][newCharger][0] + c[newCharger][curNode][0] + s[curNode];
		ChargerInPath = true;
		while (curNode != 0)
		{
			if (Direct->at(curNode).succ.size() > 0)
			{
				newNode = Direct->at(curNode).succ[0];
				if (newNode != 0)
					path.push_back(newNode);
				PathFromCharger.push_back(newNode);
				duration += c[curNode][newNode][0] + s[newNode];
				energyConsumption += c[curNode][newNode][0];
				energyConsumptionSegment += c[curNode][newNode][0];
				minCharTime = std::max(0.0, (energyConsumption - B) * r);
				curNode = newNode;

				if (duration + minCharTime > T)
				{
					if (ChargerInPath)
						AddDurationTour(&path, 0, MyCuts);
					break;

				}
				if (energyConsumptionSegment > B)
				{
					if (newNode==0)
					{
						PathFromCharger.pop_back();
						AddChargerTour(&PathFromCharger, previousCharger, 0, MyCuts);
					}
					else
					{
						AddChargerTour(&PathFromCharger, previousCharger, -1, MyCuts);
					}
					
					break;
				}
			}
			else if (Charger->at(curNode).succ.size() > 0)
			{
				newNode = Charger->at(curNode).succ[0];
				newCharger = Charger->at(curNode).succCharger[0];
				energyConsumption += c[curNode][newCharger][0] + c[newCharger][newNode][0];
				energyConsumptionSegment += c[curNode][newCharger][0];
				duration += c[curNode][newCharger][0] + c[newCharger][newNode][0];
				path.push_back(newCharger);
				if (newNode != 0)
					path.push_back(newNode);
				curNode = newNode;
				minCharTime = std::max(0.0, (energyConsumption - B) * r);
				ChargerInPath = true;
				if (duration + minCharTime > T)
				{
					AddDurationTour(&path, 0, MyCuts);
					break;
				}
				if (energyConsumptionSegment > B)
				{
					AddChargerTour(&PathFromCharger, previousCharger, newCharger, MyCuts);
					break;
				}
				PathFromCharger.clear();
				previousCharger = newCharger;
				PathFromCharger.push_back(newNode);
				energyConsumptionSegment = c[newCharger][newNode][0];
			}
			else
				break;
		}
	}

}

//ADD NO CHARGE SET
void MFGVRP_Solver::AddNoChargeSet(std::vector<int>* p, CnstrMgrPointer MyCuts) {

	int* cut = new int[p->size() + 1];
	cut[0] = 0;

	for (int i = 1; i <= p->size(); i++)
		cut[i] = p->at(i - 1);
	CMGR_AddCnstr(MyCuts, CMGR_NO_CHARGE_SET, 0, p->size(), cut, p->size() - 2);

	delete[] cut;
}

void MFGVRP_Solver::AddNoChargePath(std::vector<int>* p, CnstrMgrPointer MyCuts) {
	int* cut = new int[p->size() + 1];
	cut[0] = 0;

	for (int i = 1; i <= p->size(); i++)
		cut[i] = p->at(i - 1);

	CMGR_AddCnstr(MyCuts, CMGR_NO_CHARGE_PATH, 0, p->size(), cut, p->size() - 2);

	delete[] cut;
};


//ADD TOURNAMANT FOR CUSTOMERS
void MFGVRP_Solver::AddNoChargeTour(std::vector<int>* p, CnstrMgrPointer MyCuts) {

	int* cut = new int[p->size() + 1];
	cut[0] = 0;

	for (int i = 1; i <= p->size(); i++)
		cut[i] = p->at(i - 1);

	CMGR_AddCnstr(MyCuts, CMGR_NO_CHARGE_TOUR, 0, p->size(), cut, p->size() - 2);

	delete[] cut;
}


//ADD TOURNAMENT FOR CHARGER:
void MFGVRP_Solver::AddChargerTour(std::vector<int>* p, int in_charger, int out_charger, CnstrMgrPointer MyCuts) {

	if (in_charger >= 0 && out_charger >= 0)
	{
		int* cut = new int[p->size() + 3];
		cut[0] = 0; cut[1] = in_charger;

		for (int i = 2; i <= p->size() + 1; i++)
			cut[i] = p->at(i - 2);

		cut[p->size() + 2] = out_charger;

		CMGR_AddCnstr(MyCuts, CMGR_NO_CHARGE_DUAL_TOUR, 0, p->size() + 2, cut, p->size());

		//cut[0] = 0; cut[1] = out_charger;

		//for (int i = p->size(); i >= 2; i--)
		//	cut[i] = p->at(p->size() - i);
		//cut[p->size() + 2] = in_charger;

		//
		//CMGR_AddCnstr(MyCuts, CMGR_NO_CHARGE_DUAL_TOUR, 0, p->size() + 2, cut, p->size());


		delete[] cut;
	}
	else if (in_charger >= 0) {
		int* cut = new int[p->size() + 2];
		cut[0] = 0; cut[1] = in_charger;

		for (int i = 2; i <= p->size() + 1; i++)
			cut[i] = p->at(i - 2);

		//printf("\n");
		//for (int i = 0; i <= p->size()+1; i++)
		//{
		//	printf("%d\t", cut[i]);
		//}
		CMGR_AddCnstr(MyCuts, CMGR_NO_CHARGE_IN_TOUR, 0, p->size() + 1, cut, p->size() - 1);

		//cut[0] = 0;

		//for (int i = p->size(); i >= 1; i--)
		//	cut[i] = p->at(p->size() - i);


		//cut[p->size() + 1] = in_charger;
		////printf("\n");
		////for (int i = 0; i <= p->size() + 1; i++)
		////{
		////	printf("%d\t", cut[i]);
		////}
		//CMGR_AddCnstr(MyCuts, CMGR_NO_CHARGE_OUT_TOUR, 0, p->size() + 1, cut, p->size() - 1);
		delete[] cut;
	}
	else
	{
		int* cut = new int[p->size() + 2];
		cut[0] = 0;

		for (int i = 1; i <= p->size(); i++)
			cut[i] = p->at(i - 1);


		cut[p->size() + 1] = out_charger;

		//printf("\n");
		//for (int i = 0; i <= p->size() + 1; i++)
		//{
		//	printf("%d\t", cut[i]);
		//}
		//CMGR_AddCnstr(MyCuts, CMGR_NO_CHARGE_OUT_TOUR, 0, p->size() + 1, cut, p->size() - 1);

		//cut[0] = 0; cut[1] = out_charger;

		//for (int i = p->size()+1; i >= 2; i--)
		//	cut[i] = p->at(p->size() - i);

		//printf("\n");
		//for (int i = 0; i <= p->size() + 1; i++)
		//{
		//	printf("%d\t", cut[i]);
		//}
		//CMGR_AddCnstr(MyCuts, CMGR_NO_CHARGE_IN_TOUR, 0, p->size() + 1, cut, p->size() - 1);
		delete[] cut;


	}

}
void MFGVRP_Solver::AddFixedCharger(std::vector<int>* p, int in_charger, int out_charger, CnstrMgrPointer MyCuts) {

	if (in_charger >= 0 && out_charger >= 0)
	{
		int* cut = new int[p->size() + 3];
		cut[0] = 0; cut[1] = in_charger;

		for (int i = 2; i <= p->size() + 1; i++)
			cut[i] = p->at(i - 2);

		cut[p->size() + 2] = out_charger;

		CMGR_AddCnstr(MyCuts, CMGR_FIXED_CHARGERS, 0, p->size() + 2, cut, p->size() * 3 - 2);

		delete[] cut;
	}
	else if (in_charger >= 0) {
		int* cut = new int[p->size() + 2];
		cut[0] = 0; cut[1] = in_charger;

		for (int i = 2; i <= p->size() + 1; i++)
			cut[i] = p->at(i - 2);

		//ADD CUT
		CMGR_AddCnstr(MyCuts, CMGR_FIXED_IN_CHARGER, 0, p->size() + 1, cut, p->size() * 2 - 2);

		delete[] cut;
	}
	else
	{
		int* cut = new int[p->size() + 2];

		cut[0] = 0;

		for (int i = 1; i <= p->size(); i++)
			cut[i] = p->at(i - 1);

		cut[p->size() + 1] = out_charger;

		CMGR_AddCnstr(MyCuts, CMGR_FIXED_OUT_CHARGER, 0, p->size() + 1, cut, p->size() * 2 - 2);

		delete[] cut;


	}

}
void MFGVRP_Solver::AddFixedChargingEdge(std::vector<int>* p, int in_charger, int out_charger, CnstrMgrPointer MyCuts) {
	if (in_charger >= 0 && out_charger >= 0)
	{
		int* cut = new int[p->size() + 3];
		cut[0] = 0; cut[1] = in_charger;

		for (int i = 2; i <= p->size() + 1; i++)
			cut[i] = p->at(i - 2);

		cut[p->size() + 2] = out_charger;

		CMGR_AddCnstr(MyCuts, CMGR_FIXED_EDGES, 0, p->size() + 2, cut, p->size());

		delete[] cut;
	}
	else if (in_charger >= 0) {
		int* cut = new int[p->size() + 2];
		cut[0] = 0; cut[1] = in_charger;

		for (int i = 2; i <= p->size() + 1; i++)
			cut[i] = p->at(i - 2);

		CMGR_AddCnstr(MyCuts, CMGR_FIXED_IN_EDGE, 0, p->size() + 1, cut, p->size() - 1);

		delete[] cut;
	}
	else
	{
		int* cut = new int[p->size() + 2];

		cut[0] = 0;

		for (int i = 1; i <= p->size(); i++)
			cut[i] = p->at(i - 1);

		cut[p->size() + 1] = out_charger;

		CMGR_AddCnstr(MyCuts, CMGR_FIXED_OUT_EDGE, 0, p->size() + 1, cut, p->size() - 1);

		delete[] cut;


	}
}
void MFGVRP_Solver::AddDurationConstraint(std::vector<int>* p, int ICEV, CnstrMgrPointer MyCuts) {
	int* cut = new int[p->size() + 1];
	cut[0] = 0;

	for (int i = 1; i <= p->size(); i++)
		cut[i] = p->at(i - 1);

	switch (ICEV)
	{
	case 0:
		CMGR_AddCnstr(MyCuts, CMGR_INFEASIBLE_EV_PATH, 0, p->size(), cut, p->size() - 2);
		break;
	case 1:
		CMGR_AddCnstr(MyCuts, CMGR_INFEASIBLE_ICEV_PATH, 0, p->size(), cut, p->size() - 2);
		break;
	}


	delete[] cut;
}
void MFGVRP_Solver::AddDurationTour(std::vector<int>* p, int ICEV, CnstrMgrPointer MyCuts) {

	int* cut = new int[p->size() + 1];
	cut[0] = 0;

	for (int i = 1; i <= p->size(); i++)
		cut[i] = p->at(i - 1);

	switch (ICEV)
	{
	case 0:
		CMGR_AddCnstr(MyCuts, CMGR_DURATION_TOUR_EV, 0, p->size(), cut, p->size() - 2);
		break;
	case 1:
		CMGR_AddCnstr(MyCuts, CMGR_DURATION_TOUR_ICEV, 0, p->size(), cut, p->size() - 2);
		break;
	case 2:
		CMGR_AddCnstr(MyCuts, CMGR_DURATION_TOUR_BOTH, 0, p->size(), cut, p->size() - 2);
		break;
	}


	delete[] cut;
}
void MFGVRP_Solver::AddDurationSet(std::vector<int>* p, int ICEV, CnstrMgrPointer MyCuts) {
	int* cut = new int[p->size() + 1];
	cut[0] = 0;
	for (int i = 1; i <= p->size(); i++) {
		cut[i] = p->at(i - 1);

	}
	switch (ICEV)
	{
	case 0:
		CMGR_AddCnstr(MyCuts, CMGR_DURATION_SET_EV, 0, p->size(), cut, p->size() - 2);
		break;
	case 1:
		CMGR_AddCnstr(MyCuts, CMGR_DURATION_SET_ICEV, 0, p->size(), cut, p->size() - 2);
		break;
	case 2:
		CMGR_AddCnstr(MyCuts, CMGR_DURATION_SET_BOTH, 0, p->size(), cut, p->size() - 2);
		break;
	}


	delete[] cut;
}

void MFGVRP_Solver::NoChargePathExpr(IloExpr& expr, CnstrPointer cut) {
	expr.clear();
	int* nodes = new int[cut->IntListSize];
	for (int i = 1; i <= cut->IntListSize; i++)
		nodes[i - 1] = cut->IntList[i];

	for (int i = 0; i < cut->IntListSize - 1; i++)
		expr += x[nodes[i]][nodes[i + 1]][0] + x[nodes[i + 1]][nodes[i]][0];

	delete[] nodes;
}
void MFGVRP_Solver::NoChargeTourExpr(IloExpr& expr, CnstrPointer cut) {

	expr.clear();
	int* nodes = new int[cut->IntListSize];
	for (int i = 1; i <= cut->IntListSize; i++)
		nodes[i - 1] = cut->IntList[i];

	for (int i = 0; i < cut->IntListSize - 1; i++)
		for (int j = i + 1; j < cut->IntListSize; j++)
			expr += x[nodes[i]][nodes[j]][0];

	delete[] nodes;

}
void MFGVRP_Solver::NoChargeSetExpr(IloExpr& expr, CnstrPointer cut) {
	expr.clear();
	int* nodes = new int[cut->IntListSize];
	for (int i = 1; i <= cut->IntListSize; i++)
		nodes[i - 1] = cut->IntList[i];


	for (int i = 0; i < cut->IntListSize; i++)
		for (int j = 0; j < cut->IntListSize; j++)
			expr += x[nodes[i]][nodes[j]][0];

	delete[] nodes;
}
void MFGVRP_Solver::ChargerTourExpr(IloExpr& expr, int chargerInfo, CnstrPointer cut) {
	expr.clear();
	std::vector<int> nodes;
	std::vector<int> W;
	int r1 = 0; int r2;
	for (int i = 1; i <= cut->IntListSize; i++)
		nodes.push_back(cut->IntList[i]);




	switch (chargerInfo)
	{
	case 0:
		r1 = nodes[0];

		if (r1 == 0)
			expr += x[r1][nodes[1]][0];
		else {
			for (int i = 0; i <= n; i++)
				if (std::find(nodes.begin() + 1, nodes.end(), i) == nodes.end())
					W.push_back(i);

			for (int i = 0; i < W.size(); i++)
				for (int r = 0; r < R[W[i]][nodes[1]].size(); r++)
					if (R[W[i]][nodes[1]][r] == r1)
					{
						expr += y[W[i]][nodes[1]][r];
						break;
					}
		}
		for (int i = 1; i < nodes.size() - 1; i++)
			for (int j = i + 1; j < nodes.size(); j++)
				expr += x[nodes[i]][nodes[j]][0];
		break;
	case 1:
		r2 = nodes.back();

		for (int i = 0; i < nodes.size() - 2; i++)
			for (int j = i + 1; j < nodes.size() - 1; j++)
				expr += x[nodes[i]][nodes[j]][0];

		if (r2 == 0)
			expr += x[nodes[nodes.size() - 2]][r2][0];
		else {

			for (int i = 0; i <= n; i++)
				if (std::find(nodes.begin(), nodes.end() - 1, i) == nodes.end() - 1)
					W.push_back(i);

			for (int j = 0; j < W.size(); j++)
				for (int r = 0; r < R[nodes[nodes.size() - 2]][W[j]].size(); r++)
					if (R[nodes[nodes.size() - 2]][W[j]][r] == r2)
					{
						expr += y[nodes[nodes.size() - 2]][W[j]][r];
						break;
					}
		}
		break;
	case 2:
		r1 = nodes[0];
		r2 = nodes.back();
		for (int i = 0; i <= n; i++)
			if (std::find(nodes.begin() + 1, nodes.end() - 1, i) == nodes.end() - 1)
				W.push_back(i);

		if (nodes.size() == 2)
		{
			if (r1 == 0)
			{
				for (int i = 0; i < W.size(); i++)
				{
					auto it = std::find(R[r1][W[i]].begin(), R[r1][W[i]].end(), r2);
					if (it != R[r1][W[i]].end())
					{
						size_t index = std::distance(R[r1][W[i]].begin(), it);
						expr += y[r1][W[i]][index];
					}
				}
			}
			else
			{
				for (int i = 0; i < W.size(); i++)
				{
					auto it = std::find(R[W[i]][r2].begin(), R[W[i]][r2].end(), r1);
					if (it != R[W[i]][r2].end())
					{
						size_t index = std::distance(R[W[i]][r2].begin(), it);
						expr += y[W[i]][r2][index];
					}
				}
			}
		}
		else
		{


			if (r1 == 0)
			{
				for (int j = 1; j < nodes.size() - 1; j++)
					expr += x[r1][nodes[j]][0];
			}
			else
				for (int i = 0; i < W.size(); i++)
					for (int j = 1; j < nodes.size() - 1; j++)
						for (int r = 0; r < R[W[i]][nodes[j]].size(); r++)
							if (R[W[i]][nodes[j]][r] == r1) {
								expr += y[W[i]][nodes[j]][r];
								break;
							}

			for (int i = 1; i < nodes.size() - 2; i++)
				for (int j = i + 1; j < nodes.size() - 1; j++)
					expr += x[nodes[i]][nodes[j]][0];

			if (r2 == 0)
			{
				for (int j = 1; j < nodes.size() - 1; j++)
					expr += x[nodes[j]][r2][0];
			}
			else
				for (int i = 1; i < nodes.size() - 1; i++)
					for (int j = 0; j < W.size(); j++)
						for (int r = 0; r < R[nodes[i]][W[j]].size(); r++)
							if (R[nodes[i]][W[j]][r] == r2) {
								expr += y[nodes[i]][W[j]][r];
								break;
							}
		}


		break;
	}


}
void MFGVRP_Solver::FixedChargerExpr(IloExpr& expr, int chargerInfo, CnstrPointer cut) {
	expr.clear();
	std::vector<int> nodes;
	int r1 = 0; int r2;
	for (int i = 1; i <= cut->IntListSize; i++)
		nodes.push_back(cut->IntList[i]);




	switch (chargerInfo)
	{
	case 0:
		r1 = nodes[0];

		if (r1 == 0)
			for (int j = 1; j < nodes.size(); j++)
				expr += x[r1][nodes[j]][0];
		else {
			for (int i = 0; i <= n; i++)
				for (int j = 1; j < nodes.size(); j++)
					for (int r = 0; r < R[i][nodes[j]].size(); r++)
						if (R[i][nodes[j]][r] == r1)
						{
							expr += y[i][nodes[j]][r];
							break;
						}
		}
		for (int i = 1; i < nodes.size(); i++)
			for (int j = 1; j < nodes.size(); j++)
				expr += 2 * x[nodes[i]][nodes[j]][0];
		break;
	case 1:
		r2 = nodes.back();

		for (int i = 0; i < nodes.size() - 1; i++)
			for (int j = 0; j < nodes.size() - 1; j++)
				expr += 2 * x[nodes[i]][nodes[j]][0];

		if (r2 == 0)
			for (int i = 0; i < nodes.size() - 1; i++)
				expr += x[nodes[i]][r2][0];
		else {

			for (int i = 0; i < nodes.size() - 1; i++)
				for (int j = 0; j <= n; j++)
					for (int r = 0; r < R[nodes[i]][j].size(); r++)
						if (R[nodes[i]][j][r] == r2)
						{
							expr += y[nodes[i]][j][r];
							break;
						}
		}
		break;
	case 2:
		r1 = nodes[0];
		r2 = nodes.back();


		if (r1 == 0)
		{
			for (int j = 1; j < nodes.size() - 1; j++)
				expr += x[r1][nodes[j]][0];
		}
		else
			for (int i = 0; i <= n; i++)
				for (int j = 1; j < nodes.size() - 1; j++)
					for (int r = 0; r < R[i][nodes[j]].size(); r++)
						if (R[i][nodes[j]][r] == r1) {
							expr += y[i][nodes[j]][r];
							break;
						}

		for (int i = 1; i < nodes.size() - 1; i++)
			for (int j = 1; j < nodes.size() - 1; j++)
				expr += 3 * x[nodes[i]][nodes[j]][0];

		if (r2 == 0)
		{
			for (int j = 1; j < nodes.size() - 1; j++)
				expr += x[nodes[j]][r2][0];
		}
		else
			for (int i = 1; i < nodes.size() - 1; i++)
				for (int j = 0; j <= n; j++)
					for (int r = 0; r < R[nodes[i]][j].size(); r++)
						if (R[nodes[i]][j][r] == r2) {
							expr += y[nodes[i]][j][r];
							break;
						}
		break;
	}




}
void MFGVRP_Solver::FixedChargingEdgeExpr(IloExpr& expr, int chargerInfo, CnstrPointer cut) {

	expr.clear();
	std::vector<int> nodes;
	std::vector<int> W;
	int r1 = 0; int r2;
	for (int i = 1; i <= cut->IntListSize; i++)
		nodes.push_back(cut->IntList[i]);

	switch (chargerInfo)
	{
	case 0:
		r1 = nodes[0];
		for (int i = 0; i <= n; i++)
			if (std::find(nodes.begin() + 1, nodes.end(), i) == nodes.end())
				W.push_back(i);

		for (int i = 0; i < W.size(); i++)
		{
			if (MinCharDist[W[i]] + c[W[i]][nodes[1]][0] >= c[r1][nodes[1]][0])
				expr += x[W[i]][nodes[1]][0];
			for (int r = 0; r < R[W[i]][nodes[1]].size(); r++)
				if (c[R[W[i]][nodes[1]][r]][nodes[1]][0] >= c[r1][nodes[1]][0])
					expr += y[W[i]][nodes[1]][r];
		}

		for (int i = 1; i < nodes.size(); i++)
			for (int j = 1; j < nodes.size(); j++)
				expr += x[nodes[i]][nodes[j]][0];
		break;
	case 1:
		r2 = nodes.back();

		for (int i = 0; i <= n; i++)
			if (std::find(nodes.begin(), nodes.end()-1, i) == nodes.end()-1)
				W.push_back(i);

		for (int i = 0; i < nodes.size() - 1; i++)
			for (int j = 0; j < nodes.size() - 1; j++)
				expr += x[nodes[i]][nodes[j]][0];

		for (int j = 0; j < W.size(); j++)
		{
			if (MinCharDist[W[j]] + c[nodes[nodes.size() - 2]][W[j]][0] >= c[nodes[nodes.size() - 2]][r2][0])
				expr += x[nodes[nodes.size() - 2]][W[j]][0];
			for (int r = 0; r < R[nodes[nodes.size() - 2]][W[j]].size(); r++)
				if (c[R[nodes[nodes.size() - 2]][W[j]][r]][nodes[nodes.size() - 2]][0] >= c[r2][nodes[nodes.size() - 2]][0])
					expr += y[nodes[nodes.size() - 2]][W[j]][r];
		}
		break;
	case 2:
		r1 = nodes[0];
		r2 = nodes.back();
		for (int i = 0; i <= n; i++)
			if (std::find(nodes.begin()+1, nodes.end()-1, i) == nodes.end()-1)
				W.push_back(i);
		for (int i = 0; i < W.size(); i++)
		{
			if (MinCharDist[W[i]] + c[W[i]][nodes[1]][0] >= c[r1][nodes[1]][0])
				expr += x[W[i]][nodes[1]][0];
			for (int r = 0; r < R[W[i]][nodes[1]].size(); r++)
				if (c[R[W[i]][nodes[1]][r]][nodes[1]][0] >= c[r1][nodes[1]][0])
					expr += y[W[i]][nodes[1]][r];
		}
		r1 = r1;
		for (int i = 1; i < nodes.size() - 1; i++)
			for (int j = 1; j < nodes.size() - 1; j++)
				expr += x[nodes[i]][nodes[j]][0];
		r1 = r1;
		for (int j = 0; j < W.size(); j++)
		{
			if (MinCharDist[W[j]] + c[nodes[nodes.size() - 2]][W[j]][0] >= c[nodes[nodes.size() - 2]][r2][0])
				expr += x[nodes[nodes.size() - 2]][W[j]][0];
			for (int r = 0; r < R[nodes[nodes.size() - 2]][W[j]].size(); r++)
				if (c[R[nodes[nodes.size() - 2]][W[j]][r]][nodes[nodes.size() - 2]][0] >= c[r2][nodes[nodes.size() - 2]][0])
					expr += y[nodes[nodes.size() - 2]][W[j]][r];
		}
		r1 = r1;
		break;
	}




}
void MFGVRP_Solver::DurationTourExpr(IloExpr& expr, int ICEV, CnstrPointer cut) {

	expr.clear();
	std::vector<int> customers;
	std::vector<int> chargers;
	std::vector<int> nodes;
	std::vector<int> W;
	std::vector<std::pair<int, int>> connectedChargers;

	for (int i = 1; i <= cut->IntListSize; i++) {
		if (cut->IntList[i] <= n)
			customers.push_back(cut->IntList[i]);
		else
			chargers.push_back(cut->IntList[i]);
		nodes.push_back(cut->IntList[i]);
	}

	for (int i = 0; i <= n; i++)
		if (std::find(nodes.begin() + 1, nodes.end(), i) == nodes.end())
			W.push_back(i);


	switch (ICEV)
	{
	case 0:
		if (nodes[0] > n)
			for (int i = 0; i < W.size(); i++)
			{
				for (int j = 0; j < customers.size(); j++)
				{
					auto it = std::find(R[W[i]][customers[j]].begin(), R[W[i]][customers[j]].end(), nodes[0]);
					if (it != R[W[i]][customers[j]].end())
					{
						size_t index = std::distance(R[W[i]][customers[j]].begin(), it);
						expr += y[W[i]][customers[j]][index];
					}
				}
			}

		for (int i = 0; i < nodes.size() - 1; i++)
		{
			if (nodes[i] > n)
				continue;
			for (int j = i + 1; j < nodes.size(); j++)
			{
				if (nodes[j] > n)
					for (int t = j + 1; t < nodes.size(); t++)
					{
						if (nodes[t] > n)
							continue;
						auto it = std::find(R[nodes[i]][nodes[t]].begin(), R[nodes[i]][nodes[t]].end(), nodes[j]);
						if (it != R[nodes[i]][nodes[t]].end())
						{
							size_t index = std::distance(R[nodes[i]][nodes[t]].begin(), it);
							expr += y[nodes[i]][nodes[t]][index];
						}
					}
				else
				{
					expr += x[nodes[i]][nodes[j]][0];
					for (int r = 0; r < R[nodes[i]][nodes[j]].size(); r++)
						expr += y[nodes[i]][nodes[j]][r];
				}

			}
		}
		if (nodes.back() > n)
			for (int i = 0; i < W.size(); i++)
			{
				for (int j = 0; j < customers.size(); j++)
				{
					auto it = std::find(R[customers[j]][W[i]].begin(), R[customers[j]][W[i]].end(), nodes.back());
					if (it != R[customers[j]][W[i]].end())
					{
						size_t index = std::distance(R[customers[j]][W[i]].begin(), it);
						expr += y[customers[j]][W[i]][index];
					}
				}
			}
		break;
	case 1:
		for (int i = 0; i < customers.size() - 1; i++)
			for (int j = i + 1; j < customers.size(); j++)
			{
				expr += x[customers[i]][customers[j]][1];
			}

		break;
	case 2:
		for (int i = 0; i < customers.size() - 1; i++)
			for (int j = i + 1; j < customers.size(); j++)
			{
				expr += x[customers[i]][customers[j]][0] + x[customers[i]][customers[j]][1];
				for (int r = 0; r < R[customers[i]][customers[j]].size(); r++)
					expr += y[customers[i]][customers[j]][r];
			}

		break;
	}
}
void MFGVRP_Solver::DurationSetExpr(IloExpr& expr, int ICEV, CnstrPointer cut) {
	expr.clear();
	std::vector<int> customers;
	std::vector<int> chargers;
	std::vector<int> W;
	for (int i = 1; i <= cut->IntListSize; i++)
		if (cut->IntList[i] <= n)
			customers.push_back(cut->IntList[i]);
		else
			chargers.push_back(cut->IntList[i]);


	for (int i = 0; i <= n; i++)
		if (std::find(customers.begin() + 1, customers.end(), i) == customers.end())
			W.push_back(i);

	switch (ICEV)
	{
	case 0:


		for (int i = 0; i < customers.size(); i++)
		{
			for (int j = 0; j < customers.size(); j++)
			{
				expr += x[customers[i]][customers[j]][0];
				for (int r = 0; r < chargers.size(); r++)
				{
					auto it = std::find(R[customers[i]][customers[j]].begin(), R[customers[i]][customers[j]].end(), chargers[r]);
					if (it != R[customers[i]][customers[j]].end())
					{
						size_t index = std::distance(R[customers[i]][customers[j]].begin(), it);
						expr += 2 * y[customers[i]][customers[j]][index];
					}

				}
			}
		}

		//Add ingoing and out going through chargers:
		for (int i = 0; i < W.size(); i++)
		{
			for (int j = 0; j < customers.size(); j++)
			{
				for (int r = 0; r < chargers.size(); r++)
				{
					auto it = std::find(R[W[i]][customers[j]].begin(), R[W[i]][customers[j]].end(), chargers[r]);
					if (it != R[W[i]][customers[j]].end())
					{
						size_t index = std::distance(R[W[i]][customers[j]].begin(), it);
						expr += y[W[i]][customers[j]][index];
					}

					auto it1 = std::find(R[customers[j]][W[i]].begin(), R[customers[j]][W[i]].end(), chargers[r]);
					if (it1 != R[customers[j]][W[i]].end())
					{
						size_t index = std::distance(R[customers[j]][W[i]].begin(), it1);
						expr += y[customers[j]][W[i]][index];
					}
				}
			}
		}
		break;
	case 1:
		for (int i = 0; i < customers.size(); i++)
			for (int j = 0; j < customers.size(); j++)
				expr += x[customers[i]][customers[j]][1];


		break;
	case 2:
		for (int i = 0; i < customers.size(); i++)
			for (int j = 0; j < customers.size(); j++)
			{
				expr += x[customers[i]][customers[j]][0] + x[customers[i]][customers[j]][1];
				for (int r = 0; r < R[customers[i]][customers[j]].size(); r++)
					expr += y[customers[i]][customers[j]][r];
			}

		break;
	}

}
void MFGVRP_Solver::MTSTARExpr(IloExpr& expr, CnstrPointer cut) {



	std::vector<int> NList;
	std::vector<int> W;
	int* TList = new int[n + 1];
	int* CList = new int[n + 1];


	/* Nucleus: */
	for (int j = 1; j <= cut->IntListSize; j++)
		NList.push_back(cut->IntList[j] == n + 1 ? 0 : cut->IntList[j]);
	/* Satellites: */
	for (int j = 1; j <= cut->ExtListSize; j++)
		if (cut->ExtList[j] == n + 1)
			TList[j] = 0;
		else
			TList[j] = cut->ExtList[j];
	/* Connectors: */
	for (int j = 1; j <= cut->CListSize; j++)
		if (cut->CList[j] == n + 1)
			CList[j] = 0;
		else
			CList[j] = cut->CList[j];

	for (int i = 0; i <= n; i++)
		if (std::find(NList.begin(), NList.end(), i) == NList.end())
			W.push_back(i);

	/* Coefficients of the cut: */
	int intA = cut->A;
	int intB = cut->B;
	int intL = cut->L;
	/* Lambda=L/B, Sigma=A/B */

	/*Add the cut to the LP*/
	expr.clear();
	for (int i = 0; i < W.size(); i++)
	{
		for (int j = 1; j <= cut->IntListSize; j++)
		{
			expr += intB * (x[W[i]][cut->IntList[j]][0] + x[W[i]][cut->IntList[j]][1] + x[cut->IntList[j]][W[i]][0] + x[cut->IntList[j]][W[i]][1]);
			for (int r = 0; r < R[W[i]][cut->IntList[j]].size(); r++) expr += intB * (y[W[i]][cut->IntList[j]][r]);
			for (int r = 0; r < R[cut->IntList[j]][W[i]].size(); r++) expr += intB * (y[cut->IntList[j]][W[i]][r]);
		}
	}
	for (int i = 1; i <= cut->CListSize; i++)
	{
		for (int j = 1; j <= cut->ExtListSize; j++)
		{
			expr += -intA * (x[CList[i]][TList[j]][0] + x[TList[j]][CList[i]][0] + x[CList[i]][TList[j]][1] + x[TList[j]][CList[i]][1]);
			for (int r = 0; r < R[CList[i]][TList[j]].size(); r++) expr += -intA * (y[CList[i]][TList[j]][r] + y[TList[j]][CList[i]][r]);
		}
	}

	delete[] TList;
	delete[] CList;
}
void MFGVRP_Solver::CAPExpr(IloExpr& expr, CnstrPointer cut) {

	std::vector<int> nodes;
	expr.clear();


	for (int i = 1; i <= cut->IntListSize; ++i)
		nodes.push_back(cut->IntList[i] == n + 1 ? 0 : cut->IntList[i]);


	for (int i = 0; i < nodes.size(); i++)
		for (int j = 0; j < nodes.size(); j++)
		{
			expr += x[nodes[i]][nodes[j]][0] + x[nodes[i]][nodes[j]][1];
			for (int r = 0; r < R[nodes[i]][nodes[j]].size(); r++)
				expr += y[nodes[i]][nodes[j]][r];
		}


}
void MFGVRP_Solver::FCIExpr(IloExpr& expr, CnstrPointer cut) {
	//For FCI:
	int* Label = new int[n + 1];
	int MaxIdx;
	int MinIdx;
	int k;
	expr.clear();

	for (int j = 0; j <= n; j++) Label[j] = 0;

	MaxIdx = 0;
	for (int SubsetNr = 1; SubsetNr <= cut->ExtListSize; SubsetNr++) {

		MinIdx = MaxIdx + 1;
		MaxIdx = MinIdx + cut->ExtList[SubsetNr] - 1;

		for (int j = MinIdx; j <= MaxIdx; j++)
			Label[cut->IntList[j]] = SubsetNr;
	}

	printf("\nFCI Label\n");
	for (int i = 0; i <= n; i++)
		printf("%d\t", i);
	printf("\n");
	for (int i = 0; i <= n; i++)
		printf("%d\t", Label[i]);

	for (int i = 0; i <= n; i++)
	{
		for (int j = 0; j <= n; j++)
		{
			if (Label[i] == 0 && Label[j] != 0)
			{
				expr += x[i][j][0] + x[i][j][1];
				for (int r = 0; r < R[i][j].size(); r++) expr += y[i][j][r];
			}
			else if (Label[i] != 0 && Label[j] == 0) {
				expr += x[i][j][0] + x[i][j][1];
				for (int r = 0; r < R[i][j].size(); r++) expr += y[i][j][r];
			}
		}
	}
	for (int p = 1; p <= cut->ExtListSize; p++)
	{
		for (int i = 0; i <= n; i++)
		{
			for (int j = 0; j <= n; j++)
			{
				if (Label[i] != p && Label[j] == p)
				{
					expr += x[i][j][0] + x[i][j][1];
					for (int r = 0; r < R[i][j].size(); r++) expr += y[i][j][r];
				}
				else if (Label[i] == p && Label[j] != p)
				{
					expr += x[i][j][0] + x[i][j][1];
					for (int r = 0; r < R[i][j].size(); r++) expr += y[i][j][r];
				}
			}
		}
	}
	delete[] Label;
}
void MFGVRP_Solver::COMBExpr(IloExpr& expr, CnstrPointer cut) {
	// Book memory:
	int NoOfTeeth = cut->Key;
	int j = 0;
	int MaxIdx;
	int MinIdx;
	char** InTooth = new char* [n + 2];
	for (int i = 0; i <= n + 1; i++) InTooth[i] = new char[NoOfTeeth + 1];

	for (int Node = 0; Node <= n + 1; Node++)
		for (int Tooth = 0; Tooth <= NoOfTeeth; Tooth++)
			InTooth[Node][Tooth] = 0;


	// Read cut:
	for (int k = 1; k <= cut->IntListSize; k++)
	{
		j = cut->IntList[k];
		InTooth[j][0] = 1; /* Node j is in the handle */
	}
	for (int t = 1; t <= NoOfTeeth; t++)
	{
		//printf("\n");
		MinIdx = cut->ExtList[t];
		if (t == NoOfTeeth)
			MaxIdx = cut->ExtListSize;
		else
			MaxIdx = cut->ExtList[t + 1] - 1;
		for (int k = MinIdx; k <= MaxIdx; k++)
		{
			j = cut->ExtList[k];
			//if (j == n+1) j = 0;
			InTooth[j][t] = 1; /* Node j is in tooth t */
			//printf("%d", InTooth[j][t]);
		}
	}

	//Add cut:
	expr.clear();
	for (int i = 1; i <= n + 1; i++)
	{
		for (int j = 1; j <= n + 1; j++)
		{
			if (InTooth[i][0] != InTooth[j][0]) {
				if (i != n + 1 && j != n + 1)
				{
					expr += x[i][j][0] + x[i][j][1];
					for (int r = 0; r < R[i][j].size(); r++) expr += y[i][j][r];
				}
				else if (i == n + 1 && j < n + 1)
				{
					expr += x[0][j][0] + x[0][j][1];
					for (int r = 0; r < R[0][j].size(); r++) expr += y[0][j][r];
				}
				else if (i < n + 1 && j == n + 1)
				{
					expr += x[i][0][0] + x[i][0][1];
					for (int r = 0; r < R[i][0].size(); r++) expr += y[i][0][r];
				}

			}
			for (int t = 1; t <= NoOfTeeth; t++)
			{
				if (InTooth[i][t] != InTooth[j][t])


					if (i != n + 1 && j != n + 1)
					{
						expr += x[i][j][0] + x[i][j][1];
						for (int r = 0; r < R[i][j].size(); r++) expr += y[i][j][r];
					}
					else if (i == n + 1 && j < n + 1)
					{
						expr += x[0][j][0] + x[0][j][1];
						for (int r = 0; r < R[0][j].size(); r++) expr += y[0][j][r];
					}
					else if (i < n + 1 && j == n + 1)
					{
						expr += x[i][0][0] + x[i][0][1];
						for (int r = 0; r < R[i][0].size(); r++) expr += y[i][0][r];
					}
			}
		}

	}

	for (int Node = 0; Node <= n + 1; Node++)
		delete[] InTooth[Node];
	delete[] InTooth;
}


void MFGVRP_Solver::ConnectedComponents(UsedEdges* edges, CnstrMgrPointer MyCuts) {
	std::vector<bool> visited(edges->DirectConnections.size(), false);
	std::vector<int> component;
	//Check duration EV without chargers:
	for (int i = 1; i < edges->DirectConnections.size(); ++i) {
		for (int j = 1; j < visited.size(); j++)
			visited[j] = false;
		visited[0] = true;

		DFSUtill(i, visited, component, 0, 0, 0, s[i], false, false, &edges->Consolidated, MyCuts);
		component.clear();
	}
	//Check duration ICEV:
	if (MyCuts->Size==0)
		for (int i = 1; i < edges->ICEVConnections.size(); ++i) {
			for (int j = 1; j < visited.size(); j++)
				visited[j] = false;
			visited[0] = true;

			DFSUtill(i, visited, component, 0, 0, 0, s[i], false, true, &edges->ICEVConnections, MyCuts);
			component.clear();
		}




}
void MFGVRP_Solver::DFSUtill(int v, std::vector<bool>& visited, std::vector<int>& component, double totalWeight, double consumEnergy, double travelTime, double serviceTime, bool energy, bool ICEV, std::vector<connections>* graph, CnstrMgrPointer MyCuts) {

	visited[v] = true;

	//Sum connections
	for (int i = 0; i < component.size(); i++)
	{
		auto it = std::find(graph->at(v).pre.begin(), graph->at(v).pre.end(), component[i]);
		if (it != graph->at(v).pre.end())
		{
			size_t index = std::distance(graph->at(v).pre.begin(), it);
			totalWeight += graph->at(v).preVals[index];
		}
		it = std::find(graph->at(v).succ.begin(), graph->at(v).succ.end(), component[i]);
		if (it != graph->at(v).succ.end())
		{
			size_t index = std::distance(graph->at(v).succ.begin(), it);
			totalWeight += graph->at(v).succVals[index];
		}
	}

	component.push_back(v);


	if (energy) {
		if (MinCharDist[component[0]] + consumEnergy + MinCharDist[component.back()] > B && component.size() - 1.999 < totalWeight)
			if (TSPConcorde(&component) > B) {
				AddNoChargeSet(&component, MyCuts);
				return;
			}
	}

	if (ICEV) {
		if (c[component[0]][0][1] + travelTime + serviceTime + c[component.back()][0][1] > T && component.size() - 1.999 < totalWeight)
			if (TSPConcorde(&component,1, true) + serviceTime > T) {
				AddDurationSet(&component, 1, MyCuts);
				return;
			}
	}
	else {
		if (c[component[0]][0][0] + travelTime + serviceTime + c[component.back()][0][0] > T && component.size() - 1.999 < totalWeight) {
			double TSPval = TSPConcorde(&component, 0, true);
			if (TSPval + serviceTime > T) {
				AddDurationSet(&component, 2, MyCuts);
				return;
			}
			if (TSPval + serviceTime + std::max(0.0, (TSPval - B) * r) > T) {
				AddDurationSet(&component, 0, MyCuts);
				return;
			}
		}

	}
	//Find new connections
	for (int i = 0; i < graph->at(v).succ.size(); ++i) {
		int neighbor = graph->at(v).succ[i];
		double weight = graph->at(v).succVals[i];
		if (!visited[neighbor]) {
			consumEnergy += c[component.back()][v][0];
			if (ICEV)
				travelTime += c[component.back()][v][1];
			else 
				travelTime += c[component.back()][v][0];
			
			serviceTime += s[v];
			DFSUtill(neighbor, visited, component, totalWeight, consumEnergy, travelTime, serviceTime, energy, ICEV, graph, MyCuts);
		}
	}

	for (int i = 0; i < graph->at(v).pre.size(); ++i) {
		int neighbor = graph->at(v).pre[i];
		double weight = graph->at(v).preVals[i];
		if (!visited[neighbor]) {
			consumEnergy += c[component.back()][v][0];
			if (ICEV)
				travelTime += c[component.back()][v][1];
			else
				travelTime += c[component.back()][v][0];
			serviceTime += s[v];
			DFSUtill(neighbor, visited, component, totalWeight, consumEnergy, travelTime, serviceTime, energy, ICEV, graph, MyCuts);
		}
	}

}

