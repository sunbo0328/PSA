#ifndef _CIeeeSCCh_
#define _CIeeeSCCh_

#include "CIeeePwr.h"

typedef struct SCC_BUS {//短路电流母线计算结果
	string name;	/* 中文描述 */
	int ibs;	/* 母线号 */
	int iisland;	/* 所在岛号 */
	float kvvl;	/* 电压等级 */
	int flag;	/* 类型标志 */
	int qv;	/* 越限标志 */
	MYCOMPLEX v1;//正序电压
	MYCOMPLEX v2;//负序电压
	MYCOMPLEX v0;//零序电压
	MYCOMPLEX v_a;//A相电压
	MYCOMPLEX v_b;//B相电压
	MYCOMPLEX v_c;//C相电压
	MYCOMPLEX v;//短路线电压
	MYCOMPLEX v_pre;//故障前线电压
	float va;	/* A相幅值 */
	float vdgra;	/* A相相角 */
	float vb;	/* B相幅值 */
	float vdgrb;	/* B相相角 */
	float vc;	/* C相幅值 */
	float vdgrc;	/* C相相角 */
	float vthree;	/* 三相短路线电压幅值 */
	float athree;   /* 三相短路线电压相角 */
	int bsindex;//bus索引
}SCC_BUS;
typedef struct SCC_BRANCH {//短路电流支路计算结果
	string name;	/* 中文描述 */
	int iisland;	/* 所在岛号 */
	int ibranch;//对应IeeeModel::branch索引
	int ibs;	/* 母线号 */
	float p_i;	/* 初始有功 */
	float q_i;	/* 初始无功 */
	float p_j;	/* 初始有功 */
	float q_j;	/* 初始无功 */
	//i侧
	float i_i;	/* 三相短路电流 */
	MYCOMPLEX i3_i;/*三相短路电流*/
	MYCOMPLEX i1_i;	/* 正序电流 */
	MYCOMPLEX i2_i;	/* 负序电流 */
	MYCOMPLEX i0_i;	/* 零序电流 */
	MYCOMPLEX i_a_i;	/* A相电流 */
	MYCOMPLEX i_b_i;	/* B相电流 */
	MYCOMPLEX i_c_i;	/*C相电流 */
	//j侧
	float i_j;	/* 三相短路电流 */
	MYCOMPLEX i3_j;/*三相短路电流*/
	MYCOMPLEX i1_j;	/* 正序电流 */
	MYCOMPLEX i2_j;	/* 负序电流 */
	MYCOMPLEX i0_j;	/* 零序电流 */
	MYCOMPLEX i_a_j;	/* A相电流 */
	MYCOMPLEX i_b_j;	/* B相电流 */
	MYCOMPLEX i_c_j;	/*C相电流 */
}SCC_BRANCH;

typedef enum eABCPhase {
	eAPhase = 1,//A相
	eBPhase = 2,//A相
	eCPhase = 3//A相
}eABCPhase;
typedef class UnsymmetricalShortcircuit_OUT {//不对称故障输出
public:
	eABCPhase SpecialPhase;//特殊相。单相接地短路、一相断线时的故障相或两相相间短路、两相接地短路、两相断线故障时的正常相
	MYCOMPLEX v_fp;//故障点电压向量
	MYCOMPLEX ia_fp;//故障点A相电流
	MYCOMPLEX ib_fp;//故障点B相电流
	MYCOMPLEX ic_fp;//故障点C相电流
	MYCOMPLEX i1_fp;//故障点正序电流
	MYCOMPLEX i2_fp;//故障点负序电流
	MYCOMPLEX i0_fp;//故障点零序电流
	MYCOMPLEX va_fp;//故障点A相电压
	MYCOMPLEX vb_fp;//故障点B相电压
	MYCOMPLEX vc_fp;//故障点C相电压
	MYCOMPLEX v1_fp;//故障点正序电压
	MYCOMPLEX v2_fp;//故障点负序电压
	MYCOMPLEX v0_fp;//故障点零序电压
	vector<SCC_BUS> VBus_scc;//各节点短路后电压
	vector<SCC_BRANCH> IBranch_scc;//各支路短路电流
	void clear() {
		MYCOMPLEX temp;
		temp.r = 0;
		temp.i = 0;
		ia_fp = temp;
		ib_fp = temp;
		ic_fp = temp;
		i1_fp = temp;
		i2_fp = temp;
		i0_fp = temp;
		va_fp = temp;
		vb_fp = temp;
		vc_fp = temp;
		v1_fp = temp;
		v2_fp = temp;
		v0_fp = temp;
		VBus_scc.clear();
		IBranch_scc.clear();
	}
	UnsymmetricalShortcircuit_OUT() {
		clear();
	}
}UnsymmetricalShortcircuit_OUT;
typedef class ThreePhaseShortcircuit_OUT{//三相短路故障
public:
	MYCOMPLEX v_fp;//故障点电压向量
	MYCOMPLEX if3;//短路电流向量
	float i_fp;//故障点短路电流
	float mva_fp;//故障点短路容量
	vector<SCC_BUS> VBus_scc;//各节点短路后电压
	vector<SCC_BRANCH> IBranch_scc;//各支路短路电流
	void clear() {
		v_fp.r = 0;
		v_fp.i = 0;
		if3.r = 0;
		if3.i = 0;
		i_fp = 0;
		mva_fp = 0;
		VBus_scc.clear();
		IBranch_scc.clear();
	}
	ThreePhaseShortcircuit_OUT() {
		clear();
	}
}ThreePhaseShortcircuit_OUT;

typedef enum eUnsymmetricalShortcircuitType {//不对称故障类型
	eThreePhaseShortcircuit = 1,//三相短路
	eSinglePhaseGroundShortCircuit = 2,//单相接地短路
	eTwoPhaseShortCircuit = 3,//两相相间短路
	eTwoPhaseGroundShortCircuit = 4,//两相接地短路
	eSinglePhaseBreak = 5,//一相断线
	eTwoPhaseBreak = 6//两相断线
}eUnsymmetricalShortcircuitType;
typedef class CIeeeSCCBase :virtual public IeeeToolBase {//当前短路电流计算方法为基于潮流计算，初始化输入参数MSrc应该为状态估计或潮流计算后的计算模型
private:
	bool flag2_0;//是否计算负、零两序网络节点导纳矩阵，true则计算，false则仅计算正序
	SparseMatrixByList_GB gb1;//正序节点导纳矩阵
	//SparseMatrixByList_GB gb2;//负序节点导纳矩阵。被封上了，因为多数设备正负序阻抗相等（机组和负荷单独处理），可以用正序gb1处理后，生成GB2
	SparseMatrixByList_GB gb2;	//测试，解除注释
	SparseMatrixByList_GB gb0;//零序节点导纳矩阵
protected:
	vector<bool> nullnode0;//下表为内部节点号，true表示零序网络空节点标志（不连接任何支路零序等值电路的节点），false则为普通节点 
	//int RXInit(eUnsymmetricalShortcircuitType scc_type, int ibs_index, int jbs_index = -1);//生成节点阻抗矩阵，scc_type为不对称故障类型，ibs_index、jbs_index为bus表索引，ibs_index<0时为短路扫描；ibs_index>0时则计算指定计算母线短路电流；当为断线故障时，jbs_index必须赋有效值。
	int ThreePhaseShortcircuit(int bs_index, int nthread = 1, double rf = 0, double xf = 0);//三相短路，bs_index为bus表索引，nthread设置并行线程数，经阻抗rf、xf短路
	int SinglePhaseGroundShortCircuit(int bs_index, int nthread = 1, eABCPhase sphase = eAPhase,double rf = 0, double xf = 0);//单相接地短路，bs_index为bus表索引
	int TwoPhaseShortCircuit(int bs_index, eABCPhase sphase, double rf = 0, double xf = 0);//两相相间短路，bs_index为bus表索引
	int TwoPhaseGroundShortCircuit(int bs_index, eABCPhase sphase, double rf = 0, double xf = 0);//两相接地短路，bs_index为bus表索引
	int SinglePhaseBreak(int ibs_index,int jbs_index, eABCPhase sphase);//一相断线，ibs_index、jbs_index为bus表索引
	int TwoPhaseBreak(int ibs_index, int jbs_index, eABCPhase sphase);//两相断线，ibs_index、jbs_index为bus表索引
public:
	vector<ThreePhaseShortcircuit_OUT> scc3_out;//三相短路结果
	vector<UnsymmetricalShortcircuit_OUT> scc1_out;//单相接地短路结果
	UnsymmetricalShortcircuit_OUT scc2;//两相相间短路结果
	UnsymmetricalShortcircuit_OUT scc2_g;//两相接地短路结果
	UnsymmetricalShortcircuit_OUT m_scc1Break_out;//一相断线结果
	UnsymmetricalShortcircuit_OUT m_scc2Break_out;//两相断线结果

	string shortCircuitType;
	string sSCCType(eUnsymmetricalShortcircuitType type) {
		string stype;
		if (type == eThreePhaseShortcircuit) {
			stype = "三相短路";
			shortCircuitType = "ThreePhaseShortcircuit";
		}
		else if (type == eSinglePhaseGroundShortCircuit) {
			stype = "单相接地短路";
			shortCircuitType = "SinglePhaseGroundShortCircuit";
		}
		else if (type == eTwoPhaseShortCircuit) {
			stype = "两相相间短路";
			shortCircuitType = "TwoPhaseShortCircuit";
		}
		else if (type == eTwoPhaseGroundShortCircuit) {
			stype = "两相接地短路";
			shortCircuitType = "TwoPhaseGroundShortCircuit";
		}
		else if (type == eSinglePhaseBreak) {
			stype = "一相断线";
			shortCircuitType = "SinglePhaseBreak";
		}
		else if (type == eTwoPhaseBreak) {
			stype = "两相断线";
			shortCircuitType = "TwoPhaseBreak";
		}
		return stype;
	}
	CSparseMatrix3_KLU GB1;//正序节点导纳矩阵
	CSparseMatrix3_KLU GB2;//正序节点导纳矩阵
	CSparseMatrix3_KLU GB0;//正序节点导纳矩阵
	int GetSCCGB(int GBType);//获得短路电流计算所需正、负、零序网络节点导纳矩阵；GBType=1，则仅计算正序网络节点导纳矩阵，为其他值则同时计算正、负、零序节点导纳矩阵；先生成链表形式矩阵，再转为KLU形式矩阵
	int InitSCCIeeeModel(IeeeToolBase *MSrc);//当前短路电流计算方法为基于潮流计算，初始化输入参数MSrc应该为状态估计或潮流计算后的计算模型
	
	int SCCCalc(eUnsymmetricalShortcircuitType scc_type, int ibs_index, int jbs_index = -1);//短路电流计算，scc_type为不对称故障类型，ibs_index、jbs_index为bus表索引，ibs_index<0时为短路扫描；ibs_index>0时则计算指定计算母线短路电流；当为断线故障时，jbs_index必须赋有效值。
	int printGB(int gb1_flag = 0, int gb2_flag = 0, int gb0_flag = 0);	//主要用于测试，打印节点导纳矩阵，链表形式上三角矩阵和klu形式全矩阵，分别生成文件
	CIeeeSCCBase() {
		flag2_0 = false;
	}
}CIeeeSCCBase;



//负负荷当负荷处理，直流换流变当负荷处理，找不到次暂态电抗的机组当负荷处理？？？
typedef class CIeeeSCC :virtual public IeeeToolBase {
private:
protected:
	//int nbus;
	//int nbranch;
	//IeeeToolByList *ModelTools;
	double *NodeRX;
	SparseMatrixByList_GB S3GB;//节点导纳
	CSparseMatrix3_KLU S3GB_KLU;//节点导纳
public:
	double *If3;
	double **Id3;
	int GetNodeS3GB();//计算母线三项短路节点导纳矩阵
	int SCCModelInit(IeeeToolBase *MSrc);
	int ScanBSSCC_3(
		int ibs//ibs<0扫描所有母线三项短路，nbubs>ibs>=0扫描指定计算母线三项短路，ibs为bus逻辑号
		);

	CIeeeSCC();
	~CIeeeSCC();
	friend class ShortCircuitRatioCalc;
}CIeeeSCC;

typedef struct SetMiescr {
	string miescrname;//直流系统名称
	string stname;//受端换流站名称
	float kv;//电压等级
	float mva;//额定容量
}SetMiescr;
typedef class ShortCircuitRatioCalc  {//短路比计算类
private:
	size_t nSCRBus;//求解短路比的节点数
	string *DCCNVName;//直流名称
	int *iNBBS;//内部节点号集合
	int *iSCRBus;//外部节点号集合，即bus下标
	bool *IsDCCNV;//节点换流器是否投运
	double **iRX;//iRX[i]为i列节点阻抗矩阵元素
	double *IF3;//节点短路电流
	MYCOMPLEX **RX;//计算短路比节点的阻抗
	double **MIIF;//影响因子
	double *QF;//无功补偿
	double *MIESCR;//短路比100*v/z
	double *MIESCR2;//短路比v*v/z
	double *fscc;//短路容量100*v/z
	double *fscc2;//短路容量v*v/z
	int lv_miescr;//直流系统个数
	double *dccnvp;//节点换流变功率之和
	double *dccnvq;//节点换流变功率之和
protected:
	CIeeeSCC scc;
public:
	int GetSCRNum() {
		return lv_miescr;
	}
	vector<SetMiescr> vSetMiescr;
	int ReadSetMiescr(string fname);
	void SCRClear();
	void SCRClearAll() {
		SCRClear();
		scc.ITBClear();
	}
	int SCRInit(IeeeToolBase *MSrc);
	//int ShortCircuitRatio(const double *Bufblb1);
	int ShortCircuitRatio();
	int PrintMIESCR();

	ShortCircuitRatioCalc() {
		nSCRBus = 0;//求解短路比的节点数
		DCCNVName = NULL;//直流名称
		iNBBS = NULL;//内部节点号集合
		iSCRBus = NULL;//外部节点号集合，即bus下标
		IsDCCNV = NULL;//节点换流器是否投运
		iRX = NULL;//iRX[i]为i列节点阻抗矩阵元素
		IF3 = NULL;//节点短路电流
		RX = NULL;//计算短路比节点的阻抗
		MIIF = NULL;//影响因子
		QF = NULL;//无功补偿
		MIESCR = NULL;//短路比100*v/z
		MIESCR2 = NULL;//短路比v*v/z
		fscc = NULL;//短路容量100*v/z
		fscc2 = NULL;//短路容量v*v/z
		dccnvp = NULL;
		dccnvq = NULL;
	}
	~ShortCircuitRatioCalc() {
		SCRClear();
	}
}ShortCircuitRatioCalc;

#endif
