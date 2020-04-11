#ifndef _CIeeeSCCh_
#define _CIeeeSCCh_

#include "CIeeePwr.h"

typedef struct SCC_BUS {//��·����ĸ�߼�����
	string name;	/* �������� */
	int ibs;	/* ĸ�ߺ� */
	int iisland;	/* ���ڵ��� */
	float kvvl;	/* ��ѹ�ȼ� */
	int flag;	/* ���ͱ�־ */
	int qv;	/* Խ�ޱ�־ */
	MYCOMPLEX v1;//�����ѹ
	MYCOMPLEX v2;//�����ѹ
	MYCOMPLEX v0;//�����ѹ
	MYCOMPLEX v_a;//A���ѹ
	MYCOMPLEX v_b;//B���ѹ
	MYCOMPLEX v_c;//C���ѹ
	MYCOMPLEX v;//��·�ߵ�ѹ
	MYCOMPLEX v_pre;//����ǰ�ߵ�ѹ
	float va;	/* A���ֵ */
	float vdgra;	/* A����� */
	float vb;	/* B���ֵ */
	float vdgrb;	/* B����� */
	float vc;	/* C���ֵ */
	float vdgrc;	/* C����� */
	float vthree;	/* �����·�ߵ�ѹ��ֵ */
	float athree;   /* �����·�ߵ�ѹ��� */
	int bsindex;//bus����
}SCC_BUS;
typedef struct SCC_BRANCH {//��·����֧·������
	string name;	/* �������� */
	int iisland;	/* ���ڵ��� */
	int ibranch;//��ӦIeeeModel::branch����
	int ibs;	/* ĸ�ߺ� */
	float p_i;	/* ��ʼ�й� */
	float q_i;	/* ��ʼ�޹� */
	float p_j;	/* ��ʼ�й� */
	float q_j;	/* ��ʼ�޹� */
	//i��
	float i_i;	/* �����·���� */
	MYCOMPLEX i3_i;/*�����·����*/
	MYCOMPLEX i1_i;	/* ������� */
	MYCOMPLEX i2_i;	/* ������� */
	MYCOMPLEX i0_i;	/* ������� */
	MYCOMPLEX i_a_i;	/* A����� */
	MYCOMPLEX i_b_i;	/* B����� */
	MYCOMPLEX i_c_i;	/*C����� */
	//j��
	float i_j;	/* �����·���� */
	MYCOMPLEX i3_j;/*�����·����*/
	MYCOMPLEX i1_j;	/* ������� */
	MYCOMPLEX i2_j;	/* ������� */
	MYCOMPLEX i0_j;	/* ������� */
	MYCOMPLEX i_a_j;	/* A����� */
	MYCOMPLEX i_b_j;	/* B����� */
	MYCOMPLEX i_c_j;	/*C����� */
}SCC_BRANCH;

typedef enum eABCPhase {
	eAPhase = 1,//A��
	eBPhase = 2,//A��
	eCPhase = 3//A��
}eABCPhase;
typedef class UnsymmetricalShortcircuit_OUT {//���Գƹ������
public:
	eABCPhase SpecialPhase;//�����ࡣ����ӵض�·��һ�����ʱ�Ĺ��������������·������ӵض�·��������߹���ʱ��������
	MYCOMPLEX v_fp;//���ϵ��ѹ����
	MYCOMPLEX ia_fp;//���ϵ�A�����
	MYCOMPLEX ib_fp;//���ϵ�B�����
	MYCOMPLEX ic_fp;//���ϵ�C�����
	MYCOMPLEX i1_fp;//���ϵ��������
	MYCOMPLEX i2_fp;//���ϵ㸺�����
	MYCOMPLEX i0_fp;//���ϵ��������
	MYCOMPLEX va_fp;//���ϵ�A���ѹ
	MYCOMPLEX vb_fp;//���ϵ�B���ѹ
	MYCOMPLEX vc_fp;//���ϵ�C���ѹ
	MYCOMPLEX v1_fp;//���ϵ������ѹ
	MYCOMPLEX v2_fp;//���ϵ㸺���ѹ
	MYCOMPLEX v0_fp;//���ϵ������ѹ
	vector<SCC_BUS> VBus_scc;//���ڵ��·���ѹ
	vector<SCC_BRANCH> IBranch_scc;//��֧·��·����
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
typedef class ThreePhaseShortcircuit_OUT{//�����·����
public:
	MYCOMPLEX v_fp;//���ϵ��ѹ����
	MYCOMPLEX if3;//��·��������
	float i_fp;//���ϵ��·����
	float mva_fp;//���ϵ��·����
	vector<SCC_BUS> VBus_scc;//���ڵ��·���ѹ
	vector<SCC_BRANCH> IBranch_scc;//��֧·��·����
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

typedef enum eUnsymmetricalShortcircuitType {//���Գƹ�������
	eThreePhaseShortcircuit = 1,//�����·
	eSinglePhaseGroundShortCircuit = 2,//����ӵض�·
	eTwoPhaseShortCircuit = 3,//��������·
	eTwoPhaseGroundShortCircuit = 4,//����ӵض�·
	eSinglePhaseBreak = 5,//һ�����
	eTwoPhaseBreak = 6//�������
}eUnsymmetricalShortcircuitType;
typedef class CIeeeSCCBase :virtual public IeeeToolBase {//��ǰ��·�������㷽��Ϊ���ڳ������㣬��ʼ���������MSrcӦ��Ϊ״̬���ƻ��������ļ���ģ��
private:
	bool flag2_0;//�Ƿ���㸺������������ڵ㵼�ɾ���true����㣬false�����������
	SparseMatrixByList_GB gb1;//����ڵ㵼�ɾ���
	//SparseMatrixByList_GB gb2;//����ڵ㵼�ɾ��󡣱������ˣ���Ϊ�����豸�������迹��ȣ�����͸��ɵ�������������������gb1���������GB2
	SparseMatrixByList_GB gb2;	//���ԣ����ע��
	SparseMatrixByList_GB gb0;//����ڵ㵼�ɾ���
protected:
	vector<bool> nullnode0;//�±�Ϊ�ڲ��ڵ�ţ�true��ʾ��������սڵ��־���������κ�֧·�����ֵ��·�Ľڵ㣩��false��Ϊ��ͨ�ڵ� 
	//int RXInit(eUnsymmetricalShortcircuitType scc_type, int ibs_index, int jbs_index = -1);//���ɽڵ��迹����scc_typeΪ���Գƹ������ͣ�ibs_index��jbs_indexΪbus��������ibs_index<0ʱΪ��·ɨ�裻ibs_index>0ʱ�����ָ������ĸ�߶�·��������Ϊ���߹���ʱ��jbs_index���븳��Чֵ��
	int ThreePhaseShortcircuit(int bs_index, int nthread = 1, double rf = 0, double xf = 0);//�����·��bs_indexΪbus��������nthread���ò����߳��������迹rf��xf��·
	int SinglePhaseGroundShortCircuit(int bs_index, int nthread = 1, eABCPhase sphase = eAPhase,double rf = 0, double xf = 0);//����ӵض�·��bs_indexΪbus������
	int TwoPhaseShortCircuit(int bs_index, eABCPhase sphase, double rf = 0, double xf = 0);//��������·��bs_indexΪbus������
	int TwoPhaseGroundShortCircuit(int bs_index, eABCPhase sphase, double rf = 0, double xf = 0);//����ӵض�·��bs_indexΪbus������
	int SinglePhaseBreak(int ibs_index,int jbs_index, eABCPhase sphase);//һ����ߣ�ibs_index��jbs_indexΪbus������
	int TwoPhaseBreak(int ibs_index, int jbs_index, eABCPhase sphase);//������ߣ�ibs_index��jbs_indexΪbus������
public:
	vector<ThreePhaseShortcircuit_OUT> scc3_out;//�����·���
	vector<UnsymmetricalShortcircuit_OUT> scc1_out;//����ӵض�·���
	UnsymmetricalShortcircuit_OUT scc2;//��������·���
	UnsymmetricalShortcircuit_OUT scc2_g;//����ӵض�·���
	UnsymmetricalShortcircuit_OUT m_scc1Break_out;//һ����߽��
	UnsymmetricalShortcircuit_OUT m_scc2Break_out;//������߽��

	string shortCircuitType;
	string sSCCType(eUnsymmetricalShortcircuitType type) {
		string stype;
		if (type == eThreePhaseShortcircuit) {
			stype = "�����·";
			shortCircuitType = "ThreePhaseShortcircuit";
		}
		else if (type == eSinglePhaseGroundShortCircuit) {
			stype = "����ӵض�·";
			shortCircuitType = "SinglePhaseGroundShortCircuit";
		}
		else if (type == eTwoPhaseShortCircuit) {
			stype = "��������·";
			shortCircuitType = "TwoPhaseShortCircuit";
		}
		else if (type == eTwoPhaseGroundShortCircuit) {
			stype = "����ӵض�·";
			shortCircuitType = "TwoPhaseGroundShortCircuit";
		}
		else if (type == eSinglePhaseBreak) {
			stype = "һ�����";
			shortCircuitType = "SinglePhaseBreak";
		}
		else if (type == eTwoPhaseBreak) {
			stype = "�������";
			shortCircuitType = "TwoPhaseBreak";
		}
		return stype;
	}
	CSparseMatrix3_KLU GB1;//����ڵ㵼�ɾ���
	CSparseMatrix3_KLU GB2;//����ڵ㵼�ɾ���
	CSparseMatrix3_KLU GB0;//����ڵ㵼�ɾ���
	int GetSCCGB(int GBType);//��ö�·����������������������������ڵ㵼�ɾ���GBType=1�����������������ڵ㵼�ɾ���Ϊ����ֵ��ͬʱ����������������ڵ㵼�ɾ���������������ʽ������תΪKLU��ʽ����
	int InitSCCIeeeModel(IeeeToolBase *MSrc);//��ǰ��·�������㷽��Ϊ���ڳ������㣬��ʼ���������MSrcӦ��Ϊ״̬���ƻ��������ļ���ģ��
	
	int SCCCalc(eUnsymmetricalShortcircuitType scc_type, int ibs_index, int jbs_index = -1);//��·�������㣬scc_typeΪ���Գƹ������ͣ�ibs_index��jbs_indexΪbus��������ibs_index<0ʱΪ��·ɨ�裻ibs_index>0ʱ�����ָ������ĸ�߶�·��������Ϊ���߹���ʱ��jbs_index���븳��Чֵ��
	int printGB(int gb1_flag = 0, int gb2_flag = 0, int gb0_flag = 0);	//��Ҫ���ڲ��ԣ���ӡ�ڵ㵼�ɾ���������ʽ�����Ǿ����klu��ʽȫ���󣬷ֱ������ļ�
	CIeeeSCCBase() {
		flag2_0 = false;
	}
}CIeeeSCCBase;



//�����ɵ����ɴ���ֱ�������䵱���ɴ����Ҳ�������̬�翹�Ļ��鵱���ɴ�������
typedef class CIeeeSCC :virtual public IeeeToolBase {
private:
protected:
	//int nbus;
	//int nbranch;
	//IeeeToolByList *ModelTools;
	double *NodeRX;
	SparseMatrixByList_GB S3GB;//�ڵ㵼��
	CSparseMatrix3_KLU S3GB_KLU;//�ڵ㵼��
public:
	double *If3;
	double **Id3;
	int GetNodeS3GB();//����ĸ�������·�ڵ㵼�ɾ���
	int SCCModelInit(IeeeToolBase *MSrc);
	int ScanBSSCC_3(
		int ibs//ibs<0ɨ������ĸ�������·��nbubs>ibs>=0ɨ��ָ������ĸ�������·��ibsΪbus�߼���
		);

	CIeeeSCC();
	~CIeeeSCC();
	friend class ShortCircuitRatioCalc;
}CIeeeSCC;

typedef struct SetMiescr {
	string miescrname;//ֱ��ϵͳ����
	string stname;//�ܶ˻���վ����
	float kv;//��ѹ�ȼ�
	float mva;//�����
}SetMiescr;
typedef class ShortCircuitRatioCalc  {//��·�ȼ�����
private:
	size_t nSCRBus;//����·�ȵĽڵ���
	string *DCCNVName;//ֱ������
	int *iNBBS;//�ڲ��ڵ�ż���
	int *iSCRBus;//�ⲿ�ڵ�ż��ϣ���bus�±�
	bool *IsDCCNV;//�ڵ㻻�����Ƿ�Ͷ��
	double **iRX;//iRX[i]Ϊi�нڵ��迹����Ԫ��
	double *IF3;//�ڵ��·����
	MYCOMPLEX **RX;//�����·�Ƚڵ���迹
	double **MIIF;//Ӱ������
	double *QF;//�޹�����
	double *MIESCR;//��·��100*v/z
	double *MIESCR2;//��·��v*v/z
	double *fscc;//��·����100*v/z
	double *fscc2;//��·����v*v/z
	int lv_miescr;//ֱ��ϵͳ����
	double *dccnvp;//�ڵ㻻���书��֮��
	double *dccnvq;//�ڵ㻻���书��֮��
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
		nSCRBus = 0;//����·�ȵĽڵ���
		DCCNVName = NULL;//ֱ������
		iNBBS = NULL;//�ڲ��ڵ�ż���
		iSCRBus = NULL;//�ⲿ�ڵ�ż��ϣ���bus�±�
		IsDCCNV = NULL;//�ڵ㻻�����Ƿ�Ͷ��
		iRX = NULL;//iRX[i]Ϊi�нڵ��迹����Ԫ��
		IF3 = NULL;//�ڵ��·����
		RX = NULL;//�����·�Ƚڵ���迹
		MIIF = NULL;//Ӱ������
		QF = NULL;//�޹�����
		MIESCR = NULL;//��·��100*v/z
		MIESCR2 = NULL;//��·��v*v/z
		fscc = NULL;//��·����100*v/z
		fscc2 = NULL;//��·����v*v/z
		dccnvp = NULL;
		dccnvq = NULL;
	}
	~ShortCircuitRatioCalc() {
		SCRClear();
	}
}ShortCircuitRatioCalc;

#endif
