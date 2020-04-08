#include "CIeeeStateEstimation.h"
#include "CIeeeCA.h"
#include "MeasBaseTools.h"
#include "CIeeeSCC.h"

int main(int argc, char** argv)
{
	
	
	IeeeACPwr mypwr(1000, 2000, 200, 200);//最大1000节点，9999条支路，

	if (mypwr.StandardIeee("C:\\Users\\sunbo\\Desktop\\1.PSASource\\IEEE\\005ieee.dat") < 0) {
		ErrorMessage << "读取数据文件" << "C:\\Users\\sunbo\\Desktop\\1.PSASource\\IEEE\\005ieee.dat" << "失败！" << endl;
		return -1;
	}
	mypwr.InitPwrGlogbal();
	if (mypwr.ACPwr(DCFLAG_EQ) < 0) {
		ErrorMessage << "潮流计算失败！" << endl;
	}
	if (mypwr.BSEquality(0) < 0) {
		ErrorMessage << "BSEquality计算失败！" << endl;
	}
	stringstream stemp;
	stemp.str("");
	stemp.clear();
	stemp << "C:\\Users\\sunbo\\Desktop\\1.PSASource\\debug\\" << "DPF_IEEE_bus" << mypwr.nbus << "_branch" << mypwr.nbranch;
	string fname = stemp.str();
	mypwr.printTOPO(fname);


	CIeeeSCCBase myscc;
	if (myscc.InitSCCIeeeModel(&mypwr) < 0)
	{
		ErrorMessage << "error" << endl;
	}
	//myscc.SCCCalc(eThreePhaseShortcircuit, 1);
	//myscc.SCCCalc(eSinglePhaseGroundShortCircuit, 1);
	//myscc.SCCCalc(eTwoPhaseShortCircuit, 1);
	myscc.SCCCalc(eTwoPhaseGroundShortCircuit, 1);

	return 1;
}

int main2(int argc, char** argv)
{
	APPControlParameter message;
	if (argc < 5) {
		cout << "参数错误，示例：PSACalcDemo appCPConfig.json se(dpf|dpf1) 1(0则忽略提示信息) 1(0则不回写实时库)" << endl << endl << endl;
		message.PrintJSONDemo("NULL");
		return 1;
	}
	string fname = argv[1];
	string sCType = argv[2];
	AlarmMessageDebugOutFlag = atoi(argv[3]);
	int updatefalg = atoi(argv[4]);
	//读取进程控制参数配置文件
	if (message.InitAPPCP(fname) < 0) {
		return -1;
	}
	if(message.GetMode()!="ieee"){
		ErrorMessage << "模型和数据源类型选择错误，mode应=ieee" << endl;
		return -1;
	}

	string LogDir;
	//检查并创建运行路径
	LogDir.clear();
	LogDir = LogDir + getenv("HOME") + "/var/pasdb/log/";
	if (MakeFolder(LogDir.c_str()) < 0) {
		AlarmMessage << "创建路径" << LogDir << "失败" << endl;
		return -1;
	}
	string DataFile = message.GetFData();
	CACPInitByJson cacp;
	if (cacp.InitCACPByJson(message.GetCPFName()) < 0) {
		AlarmMessage << "获取潮流计算控制参数失败！" << endl;
		return -1;
	}

	IeeeACPwr  ACDpf;
	ACDpf.init(1000, 9999, 1000, 1000);
	if (ACDpf.StandardIeee(DataFile.c_str()) < 0) {
		AlarmMessage << "读取数据文件" << DataFile.c_str() << "失败！" << endl;
		return -1;
	}
	ACDpf.InitPwrGlogbal(cacp.dpfcp);
	if (ACDpf.ACPwr(DCFLAG_EQ) < 0) {
		AlarmMessage << "潮流计算失败！" << endl;
		return -1;
	}
	if (message.GetLogFlag()) {
		string logname = LogDir + "DPF";
		if (cacp.dpfcp.pqwr == NEWTONPWR) {
			if (ACDpf.GetJocabi() < 0) {
				ErrorMessage << "计算雅可比矩阵失败！" << endl;
				return -1;
			}
			ACDpf.PrintJocbi(logname);
		}
		ACDpf.printBSEqu(logname);
		ACDpf.printPQV(logname);
		ACDpf.printTOPO(logname);
	}

	CIeeeCA ca;
	ca.InitCAGlobal(cacp.cacp, cacp.dpfcp);//ca控制参数
	if (ca.CAInit(&ACDpf, 0, NULL) < 0) {
		AlarmMessage << "静态安全分析初始化失败！" << endl;
		return -1;
	}
	if (ca.CA_AC_N_1(1) < 0) {
		ErrorMessage << "ca.SA_AC_N_1错误！" << endl;
		return -1;
	}
	
	return 1;
}


int main_test()
{
	double _H[9] = { 1,7,1,3,5,1,-4,6,1 };
	double dWl[3] = { 0,2,8 };
	double dWy[3] = { 0,0,0 };
	int dWySize = 3;
	int numE = 0;
	for (int i = 0; i < dWySize; i++)//列
	{
		numE = i* dWySize;
		for (int j = 0; j < dWySize; j++)//行
		{
			dWy[j] += _H[numE + j] * dWl[i];
		}
	}
	for (int i = 0; i < dWySize; i++)
	{
		cout << dWy[i] << endl;
	}
	system("pause");
	return -1;

	CSparseMatrix3_KLU A;
	A.CSparseMatrix3Init(eRealNum, eSaveByColumn, 16, 4, 4);
	A.Ap[0] = 0;
	A.Ap[1] = 4;
	A.Ap[2] = 8;
	A.Ap[3] = 12;
	A.Ap[4] = 16;
	for (int i = 0; i < 16; i++)
	{
		if (i < 4) {
			A.Ai[i] = i;
		}
		else if (i < 8 && i >= 4) {
			A.Ai[i] = i - 4;
		}
		else if (i < 12 && i >= 8) {
			A.Ai[i] = i - 8;
		}
		else if (i < 16 && i >= 12) {
			A.Ai[i] = i - 12;
		}
	}
	A.Ax[0] = 2;
	A.Ax[1] = 1;
	A.Ax[2] = 2;
	A.Ax[3] = 0;
	A.Ax[4] = 2;
	A.Ax[5] = 3;
	A.Ax[6] = 4;
	A.Ax[7] = 0;
	A.Ax[8] = 3;
	A.Ax[9] = 5;
	A.Ax[10] = 7;
	A.Ax[11] = 0;
	A.Ax[12] = 0;
	A.Ax[13] = 0;
	A.Ax[14] = 0;
	A.Ax[15] = 1;
	A.Inversion();
	double *_A = A.RetIMp();
	cout << setiosflags(ios::left) << setw(16) << _A[0] << setw(16) << _A[1] << setw(16) << _A[2] << setw(16) << _A[3] << endl;
	cout << setiosflags(ios::left) << setw(16) << _A[4] << setw(16) << _A[5] << setw(16) << _A[6] << setw(16) << _A[7] << endl;
	cout << setiosflags(ios::left) << setw(16) << _A[8] << setw(16) << _A[9] << setw(16) << _A[10] << setw(16) << _A[11] << endl;
	cout << setiosflags(ios::left) << setw(16) << _A[12] << setw(16) << _A[13] << setw(16) << _A[14] << setw(16) << _A[15] << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			double value = 0;
			cout << "[" << i << "][" << j << "]=";
			for (int ii = A.Ap[j]; ii < A.Ap[j + 1]; ii++)
			{
				value = value + (_A[i + 4 * A.Ai[ii]] * A.Ax[ii]);
			}
			cout << value << endl;
		}

	}
	system("pause");
	return -1;

	vector<OffEquipment> v1;
	vector<OffEquipment> v2;
	OffEquipment oetemp;
	oetemp.iisland = 999;
	oetemp.index = 100;
	oetemp.type = 333;
	strncpy(oetemp.name, "ceshi", 128);
	for (int i = 0; i < 100; i++)
	{
		v1.push_back(oetemp);
	}
	v2 = v1;
	for (int i = 0; i < v2.size(); i++)
	{
		cout << i << " " << v2[i].iisland << " " << v2[i].index << " " << v2[i].type << " " << v2[i].name << endl;
	}

#ifdef __WINDOWS__
	system("pause");
#endif

	return 1;
	int size = 99;
	int nthread = 100;
	const int nCPUCore = omp_get_num_procs();
	const int MaxThreads = 2 * nCPUCore;
	cout << "            nCPUCore[" << nCPUCore << "],MaxThreads[" << MaxThreads << "]" << endl;
	int omp_threads = nthread > MaxThreads ? MaxThreads : nthread;
	cout << "            omp_threads[" << omp_threads << "]" << endl;
	omp_lock_t mylock;
	omp_init_lock(&mylock);
#pragma omp parallel for
	for (int i = 0; i<size; i++)
	{
		double ii = i*i*0.998;
		int nthreads = omp_get_num_threads();//线程总数
		int id = omp_get_thread_num();//线程号
		omp_set_lock(&mylock);
		cout << "i=" << i << "  thread" << id << " ii=" << ii << " nthreads=" << nthreads << endl;
		omp_unset_lock(&mylock);
	}
	omp_destroy_lock(&mylock);
#ifdef __WINDOWS__
	system("pause");
#endif

	return 1;
	IeeeBranch b1[20];
	IeeeBranch b2[13];
	for (int i = 0; i < 20; i++)
	{
		sprintf(b1[i].name, "支路%d", i);
		b1[i].id = i;
		b1[i].r = i*0.001;
		b1[i].x = i*0.002;
		b1[i].amprating = 100 * i;
		cout << b1[i].id << " " << b1[i].name << " " << b1[i].r << " " << b1[i].x << " " << b1[i].amprating << endl;
	}
	cout << endl;
	vector<int> voff;
	voff.push_back(0);
	voff.push_back(1);
	voff.push_back(5);
	voff.push_back(9);
	voff.push_back(10);
	voff.push_back(11);
	voff.push_back(12);
	IeeeBranch *sdp = b2;
	IeeeBranch *ssp = b1;
	int esp = 0, num = 0;
	for (int i = 0; i < voff.size(); i++)
	{
		esp = voff[i];
		if (i == 0) {
			num = esp;
		}
		else {
			num = (voff[i] - voff[i - 1] - 1);
		}
		if (num > 0) {
			memcpy((void *)sdp, (void *)ssp, num * sizeof(IeeeBranch));
		}
		sdp = sdp + num;
		ssp = b1 + esp + 1;
	}
	num = 20 - esp - 1;
	if (num > 0) {
		memcpy((void *)sdp, (void *)ssp, num * sizeof(IeeeBranch));
	}
	for (int i = 0; i < 13; i++)
	{
		cout << b2[i].id << " " << b2[i].name << " " << b2[i].r << " " << b2[i].x << " " << b2[i].amprating << endl;
	}
#ifdef __WINDOWS__
	system("pause");
#endif
	return 1;

	MeasBaseTools measbase;
	measbase.ReadMeasBase("measbase.log");
	float base = 0, percent = 0;int  vl = 0;
	measbase.GetBasePercent(0, 535, vl, base, percent);
	cout << base << " " << percent << endl;
#ifdef __WINDOWS__
	system("pause");
#endif
	return 1;
	IeeeACWLSSE se;
	vector<int> a;
	a.resize(1000,999);
	for (size_t i = 0; i < a.size(); i++) {
		cout << a.size() << " " << a[i] << endl;
	}
#ifdef __WINDOWS__
	system("pause");
#endif
	return 1;
}
