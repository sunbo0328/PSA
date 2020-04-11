#include "CIeeeSCC.h"

int CIeeeSCCBase::GetSCCGB(int GBType)
{
	if (!IsReady()) {
		AlarmMessage << "Ieee模型未准备完毕！" << endl;
		return -1;
	}
	TIMING jsks("CIeeeSCC::GetSCCGB");
	int it = 0;
	//double r = 0 , x = 0, x_r = 0;// , bch;
	//double xfr = 0, xfx = 0;
	//double g2b2 = 0;
	double g = 0, b = 0, bch2 = 0, t = 0;
	double g0 = 0, b0 = 0, bch0_2 = 0;
	int i = 0, j = 0;
	LISTELE * iigb = NULL;
	LISTELE * jjgb = NULL;
	LISTELE * ijgb = NULL;
	AlarmMessage << "节点导纳矩阵维数(" << iNAllACNode << "*" << iNAllACNode << ")" << endl;
	if (GBType != 1) {
		flag2_0 = true;
	}
	gb1.init(iNAllACNode);
	if (flag2_0) {
		gb0.init(iNAllACNode);
	}
	nullnode0.clear();
	nullnode0.resize(iNAllACNode, false);
	
	for (it = 0; it < nbranch; it++)
	{
		if (branch[it].branchtype == BRANCHTYPE_LN) {//线路
			i = branch[it].i;
			j = branch[it].j;
			i = bus[i].ibs;//转换内部节点
			j = bus[j].ibs;//转换内部节点
			g = branch[it].g;
			b = branch[it].b;
			bch2 = branch[it].bch2;///10000;//1/2线路电纳
			//branch[it].bch2 = 0;
			//正序
			iigb = gb1.FindOrAdd(i, i);//gb1为对称矩阵，仅保存上三角
			jjgb = gb1.FindOrAdd(j, j);
			ijgb = gb1.FindOrAdd(i, j);
			//节点i,j自阻抗
			//节点i,j自阻抗
			iigb->real += g;
			iigb->imag += (b + bch2);
			jjgb->real += g;
			jjgb->imag += (b + bch2);
			//i,j互阻抗
			ijgb->real -= g;
			ijgb->imag -= b;
#ifdef _DEBUG
			if (gb1.EleIsNanOrInf(iigb) || gb1.EleIsNanOrInf(jjgb) || gb1.EleIsNanOrInf(ijgb)) {
				AlarmMessage << "ln_branch[" << it << "],name(" << branch[it].name << "),g(" << branch[it].g << "),b(" << branch[it].b << ")" << endl;
				return -1;
			}
#endif
			//零序
			if (flag2_0) {
				nullnode0[i] = true;
				nullnode0[j] = true;
				g0 = branch[it].g0;
				b0 = branch[it].b0;
				bch0_2 = branch[it].bch0_2;///10000;//1/2线路电纳
				iigb = gb0.FindOrAdd(i, i);//gb0为对称矩阵，仅保存上三角
				jjgb = gb0.FindOrAdd(j, j);
				ijgb = gb0.FindOrAdd(i, j);
				//节点i,j自阻抗
				//节点i,j自阻抗
				iigb->real += g0;
				iigb->imag += (b0 + bch0_2);
				jjgb->real += g0;
				jjgb->imag += (b0 + bch0_2);
				//i,j互阻抗
				ijgb->real -= g0;
				ijgb->imag -= b0;
#ifdef _DEBUG
				if (gb0.EleIsNanOrInf(iigb) || gb0.EleIsNanOrInf(jjgb) || gb0.EleIsNanOrInf(ijgb)) {
					AlarmMessage << "零序ln_branch[" << it << "],name(" << branch[it].name << "),g0(" << branch[it].g << "),b0(" << branch[it].b << ")" << endl;
					return -1;
				}
#endif
			}
		}
		else if (branch[it].branchtype == BRANCHTYPE_LNS) {//串补
			i = branch[it].i;
			j = branch[it].j;
			i = bus[i].ibs;//转换内部节点
			j = bus[j].ibs;//转换内部节点
			//正序
			b = branch[it].b;
			iigb = gb1.FindOrAdd(i, i);//gb1为对称矩阵，仅保存上三角
			jjgb = gb1.FindOrAdd(j, j);
			ijgb = gb1.FindOrAdd(i, j);
			//节点i,j自阻抗
			//节点i,j自阻抗
			iigb->imag += b;
			jjgb->imag += b;
			//i,j互阻抗
			ijgb->imag -= b;
#ifdef _DEBUG
			if (gb1.EleIsNanOrInf(iigb) || gb1.EleIsNanOrInf(jjgb) || gb1.EleIsNanOrInf(ijgb)) {
				AlarmMessage << "lns_branch[" << it << "],name(" << branch[it].name << "),b(" << branch[it].b << ")" << endl;
				return -1;
			}
#endif
			//零序
			if (flag2_0) {
				nullnode0[i] = true;
				nullnode0[j] = true;
				b0 = branch[it].b0;
				iigb = gb0.FindOrAdd(i, i);//gb0为对称矩阵，仅保存上三角
				jjgb = gb0.FindOrAdd(j, j);
				ijgb = gb0.FindOrAdd(i, j);
				//节点i,j自阻抗
				//节点i,j自阻抗
				iigb->imag += b0;
				jjgb->imag += b0;
				//i,j互阻抗
				ijgb->imag -= b0;
#ifdef _DEBUG
				if (gb0.EleIsNanOrInf(iigb) || gb0.EleIsNanOrInf(jjgb) || gb0.EleIsNanOrInf(ijgb)) {
					AlarmMessage << "lns_branch[" << it << "],name(" << branch[it].name << "),b(" << branch[it].b << ")" << endl;
					return -1;
				}
#endif
			}
		}
		else if (branch[it].branchtype == BRANCHTYPE_XF) {//变压器支路//t在一端，i端
			i = branch[it].i;//高压侧，变压器参数折算到高压侧，因此，i侧为标准侧，j侧为非标准差（状态估计书中的叫法）?????????
			j = branch[it].j;//低压侧
			i = bus[i].ibs;//转换内部节点
			j = bus[j].ibs;//转换内部节点
			//正序
			g = branch[it].g;
			b = branch[it].b;
			t = branch[it].t;
			if (t<0.6 || t>1.4) {
				AlarmMessage << "xf_branch[" << it << "],name(" << branch[it].name << "),g(" << g << "),b(" << b << "),t(" << t << "),变比异常！" << endl;
				t = 1;
			}
			iigb = gb1.FindOrAdd(i, i);
			jjgb = gb1.FindOrAdd(j, j);
			ijgb = gb1.FindOrAdd(i, j);
			//i端为非标准侧，j端标准侧//节点i,j自阻抗(pwr_new中也这么处理的，fd_yy函数中tj=1)
			//branch[it].t = 1;
			iigb->real += (g / (t*t));	//参数为变比为1侧的参数，须折算到高压侧
			iigb->imag += (b / (t*t));
			jjgb->real += (g);
			jjgb->imag += (b);
			//i,j互阻抗
			ijgb->real -= (g / t);
			ijgb->imag -= (b / t);
#ifdef _DEBUG
			if (gb1.EleIsNanOrInf(iigb) || gb1.EleIsNanOrInf(jjgb) || gb1.EleIsNanOrInf(ijgb)) {
				AlarmMessage << "xf_branch[" << it << "],name(" << branch[it].name << "),g(" << branch[it].g << "),b(" << branch[it].b << "),t(" << branch[it].t << ")" << endl;
				return -1;
			}
#endif
			//零序
			if (flag2_0) {
				g0 = branch[it].g0;
				b0 = branch[it].b0;
				//t = 1;// branch[it].t;//零序变压器等值电路变比为1，即没有对地并联支路
				iigb = gb0.FindOrAdd(i, i);
				jjgb = gb0.FindOrAdd(j, j);
				ijgb = gb0.FindOrAdd(i, j);
				//i端为非标准侧，j端标准侧//节点i,j自阻抗(pwr_new中也这么处理的，fd_yy函数中tj=1)
				//branch[it].t = 1;
				if (branch[it].windingcmode == XF2_DYN) {//J侧对地阻抗
					iigb->real += EPSINON;
					iigb->imag += EPSINON;
					jjgb->real += g0;
					jjgb->imag += b0;
					//i,j互阻抗
					ijgb->real -= EPSINON;
					ijgb->imag -= EPSINON;
				}
				else if (branch[it].windingcmode == XF2_YND) {//I侧对地阻抗
					iigb->real += g0;
					iigb->imag += b0;
					jjgb->real += EPSINON;
					jjgb->imag += EPSINON;
					//i,j互阻抗
					ijgb->real -= EPSINON;
					ijgb->imag -= EPSINON;
				}
				else if (branch[it].windingcmode == XF2_YNYN) {//IJ侧支路阻抗
					nullnode0[i] = true;
					nullnode0[j] = true;
					iigb->real += g0;
					iigb->imag += b0;
					jjgb->real += g0;
					jjgb->imag += b0;
					//i,j互阻抗
					ijgb->real -= g0;
					ijgb->imag -= b0;
				}
				else if (branch[it].windingcmode == XF2_DEFAULT) {//IJ侧支路阻抗无穷大，即导纳无穷小
					iigb->real += EPSINON;
					iigb->imag += EPSINON;
					jjgb->real += EPSINON;
					jjgb->imag += EPSINON;
					//i,j互阻抗
					ijgb->real -= EPSINON;
					ijgb->imag -= EPSINON;
				}
				else if (branch[it].windingcmode == XF3_D) {//J侧对地阻抗（三卷变支路j侧为虚拟节点）
					iigb->real += EPSINON;
					iigb->imag += EPSINON;
					jjgb->real += g0;
					jjgb->imag += b0;
					//i,j互阻抗
					ijgb->real -= EPSINON;
					ijgb->imag -= EPSINON;
				}
				else if (branch[it].windingcmode == XF3_Y) {//IJ侧支路阻抗无穷大，即导纳无穷小
					iigb->real += EPSINON;
					iigb->imag += EPSINON;
					jjgb->real += EPSINON;
					jjgb->imag += EPSINON;
					//i,j互阻抗
					ijgb->real -= EPSINON;
					ijgb->imag -= EPSINON;
				}
				else if (branch[it].windingcmode == XF3_YN) {//IJ侧支路阻抗
					nullnode0[i] = true;
					nullnode0[j] = true;
					iigb->real += g0;
					iigb->imag += b0;
					jjgb->real += g0;
					jjgb->imag += b0;
					//i,j互阻抗
					ijgb->real -= g0;
					ijgb->imag -= b0;
				}
				else {//类型未指定，默认//IJ侧支路阻抗
					nullnode0[i] = true;
					nullnode0[j] = true;
					iigb->real += g0;
					iigb->imag += b0;
					jjgb->real += g0;
					jjgb->imag += b0;
					//i,j互阻抗
					ijgb->real -= g0;
					ijgb->imag -= b0;
				}
#ifdef _DEBUG
				if (gb0.EleIsNanOrInf(iigb) || gb0.EleIsNanOrInf(jjgb) || gb0.EleIsNanOrInf(ijgb)) {
					AlarmMessage << "xf_branch[" << it << "],name(" << branch[it].name << "),g(" << branch[it].g << "),b(" << branch[it].b << "),t(" << branch[it].t << ")" << endl;
					return -1;
				}
#endif
			}
		}
		else {
			AlarmMessage << "支路：branch[" << it << "](" << branch[it].name << ")->branchtype=" << branch[it].branchtype << "未处理！！！" << endl;
		}
	}
	//cout<<endl<<endl;
	for (it = 0; it < nbus; it++)//并联电容器、电抗器、负荷、机组
	{
		//if (DCBus_DCIsland[it] > 0)continue;//直流母线
		if (bus[it].nodetype == NODE_DCBS)continue;//直流母线，正常是不会有并联容抗器的
		if (bus[it].nacbranch == 0)continue;//无交流支路，即该节点为换流变侧交流节点，通过换流器与主网有连接，但换流变开关断开
		i = bus[it].ibs;
		//正序
		iigb = gb1.FindOrAdd(i, i);
		//并联电容器/电抗器
		iigb->real += bus[it].blg;
		iigb->imag += bus[it].blb;	
		if (bus[it].nld > 0 || bus_ndccnv[it] > 0) {//负荷恒阻抗模型
			iigb->real += bus[it].ldg;
			iigb->imag += bus[it].ldb;
		}
		if (bus[it].nun > 0) {//机组次暂态电抗，没参数则采用恒阻抗模型			
			if (fabs(bus[it].ungd) > EPSINON) {//机组
				//AlarmMessage << "[bs|" << bus[it].id << "].ungd=" << bus[it].ungd << "!" << endl;
				iigb->imag += bus[it].ungd;
			}
			else {
				//if (fabs(bus[it].unp) > 1.0 || fabs(bus[it].unq) > 1.0)
				{
					//AlarmMessage << "匹配不上的机组按负负荷处理！" << endl;
					if (fabs(bus[i].vc) > EPSINON) {
						iigb->real += (-bus[it].unp*_wbase) / (bus[it].vc*bus[it].vc);
						iigb->imag += (bus[it].unq*_wbase) / (bus[it].vc*bus[it].vc);
					}
					else if (fabs(bus[i].v) > EPSINON) {
						iigb->real += (-bus[it].unp*_wbase) / (bus[it].v*bus[it].v);
						iigb->imag += (bus[it].unq*_wbase) / (bus[it].v*bus[it].v);
					}
					else {
						iigb->real += (-bus[it].unp*_wbase);
						iigb->imag += (bus[it].unq*_wbase);
					}
				}
			}
		}
#ifdef _DEBUG
		if (gb1.EleIsNanOrInf(iigb)) {
			return -1;
		}
#endif
		//零序
		if (flag2_0) {
			iigb = gb0.FindOrAdd(i, i);
			//并联电容器/电抗器
			iigb->real += bus[it].blg;
			iigb->imag += bus[it].blb;
		}
	}
	gb1.SetIsReady(true);
	if (GB1.InitBySCM_Column_KLU(&gb1, eComplexNum2, 0) < 0) {
		ErrorMessage << "KLU节点导纳矩阵GB1初始化错误！！！" << endl;
		return -1;
	}
	if (flag2_0) {
		gb0.SetIsReady(true);	//零序导纳矩阵已经在上面处理完成，直接构建阻抗矩阵
		if (GB0.InitBySCM_Column_KLU(&gb0, eComplexNum2, 0) < 0) {
			ErrorMessage << "KLU节点导纳矩阵GB0初始化错误！！！" << endl;
			return -1;
		}
		//修正后形成负序节点导纳矩阵，用正序导纳矩阵构建负序导纳矩阵，并求逆得负序阻抗矩阵
		for (it = 0; it < nbus; it++)//并联电容器、电抗器、负荷、机组
		{
			//if (DCBus_DCIsland[it] > 0)continue;//直流母线
			if (bus[it].nodetype == NODE_DCBS)continue;//直流母线，正常是不会有并联容抗器的
			if (bus[it].nacbranch == 0)continue;//无交流支路，即该节点为换流变侧交流节点，通过换流器与主网有连接，但换流变开关断开
			i = bus[it].ibs;
			iigb = gb1.FindOrAdd(i, i);
			if (bus[it].nld > 0 || bus_ndccnv[it] > 0) {//负荷恒阻抗模型Z2=0.35*Z1
				iigb->real -= 0.65*bus[it].ldg;
				iigb->imag -= 0.65*bus[it].ldb;
			}
			if (bus[it].nun > 0) {//机组次暂态电抗，没参数则采用恒阻抗模型
				if (fabs(bus[it].ungd2) > EPSINON) {//提供了负序机组电抗
					if (fabs(bus[it].ungd) > EPSINON) {//机组
						iigb->imag -= bus[it].ungd;
					}
					else {
						if (fabs(bus[i].vc) > EPSINON) {
							iigb->real -= (-bus[it].unp*_wbase) / (bus[it].vc*bus[it].vc);
							iigb->imag -= (bus[it].unq*_wbase) / (bus[it].vc*bus[it].vc);
						}
						else if (fabs(bus[i].v) > EPSINON) {
							iigb->real -= (-bus[it].unp*_wbase) / (bus[it].v*bus[it].v);
							iigb->imag -= (bus[it].unq*_wbase) / (bus[it].v*bus[it].v);
						}
						else {
							iigb->real -= (-bus[it].unp*_wbase);
							iigb->imag -= (bus[it].unq*_wbase);
						}
					}
					iigb->imag += bus[it].ungd2;
				}
				else {//没有提供负序机组电抗
					if (fabs(bus[it].ungd) > EPSINON) {//机组，负序电抗默认等于正序
					}
					else {
						if (fabs(bus[i].vc) > EPSINON) {
							iigb->real -= 0.65*(-bus[it].unp*_wbase) / (bus[it].vc*bus[it].vc);
							iigb->imag -= 0.65*(bus[it].unq*_wbase) / (bus[it].vc*bus[it].vc);
						}
						else if (fabs(bus[i].v) > EPSINON) {
							iigb->real -= 0.65*(-bus[it].unp*_wbase) / (bus[it].v*bus[it].v);
							iigb->imag -= 0.65*(bus[it].unq*_wbase) / (bus[it].v*bus[it].v);
						}
						else {
							iigb->real -= 0.65*(-bus[it].unp*_wbase);
							iigb->imag -= 0.65*(bus[it].unq*_wbase);
						}
					}
				}
			}

#ifdef _DEBUG
			if (gb1.EleIsNanOrInf(iigb)) {
				return -1;
			}
#endif
		}
		if (GB2.InitBySCM_Column_KLU(&gb1, eComplexNum2, 0) < 0) {
			ErrorMessage << "KLU节点导纳矩阵GB2初始化错误！！！" << endl;
			return -1;
		}
	}

	return 1;
}
int CIeeeSCCBase::InitSCCIeeeModel(IeeeToolBase *MSrc)
{
	if (IeeeToolBaseCopy(MSrc) < 0) {
		AlarmMessage << "拷贝模型IeeeToolBaseCopy错误!" << endl;
		return -1;
	}

	return 1;
}
int CIeeeSCCBase::ThreePhaseShortcircuit(int bs_index, int nthread, double rf, double xf)//单一三相短路或三相短路接地故障
{
	if (nbus == 1) {
		ErrorMessage << "就1个节点，开玩笑？" << endl;
		return -1;
	}
	int NBUS = 0;//扫描母线总数
	int sNum = 0;//扫描起始位置
	int eNum = 0;//扫描结束位置
	double *ZZ1_bs = NULL;
	const size_t nrow = GB1.GetIDim();//节点阻抗矩阵行数
	if (bs_index < 0) {//扫描所有母线三项短路
		NBUS = nbus;
		sNum = 0;
		for (int i = 0; i < NBUS; i++)	//初始化节点计算结果vector和支路计算结果vector的维数
		{
			scc3_out[i].VBus_scc.resize(NBUS);
			scc3_out[i].IBranch_scc.resize(nbranch);
		}
		if (GB1.Inversion() < 0) {
			ErrorMessage << "矩阵求逆失败！" << endl;
			return -1;
		}
	}
	else if (bs_index>nbus) {
		ErrorMessage << "指定母线下表[" << bs_index << "]无效！" << endl;
		return -1;
	}
	else {
		NBUS = 1;
		sNum = bs_index;
		scc3_out[sNum].VBus_scc.resize(nrow);	//初始化节点计算结果vector的维数
		scc3_out[sNum].IBranch_scc.resize(nbranch);	//初始化支路计算结果vector的维数
		if (GB1.Inversion(bus[bs_index].ibs, &ZZ1_bs) < 0) {
			ErrorMessage << "求解bus[" << bs_index << "][bs|" << bus[bs_index].id << "](" << bus[bs_index].name << ")自阻抗互阻抗失败!" << endl;
			return -1;
		}
	}
	eNum = sNum + NBUS;
	int *inb_iwb = GetNBI_WBI();
	//设置并行线程
	const int nCPUCore = omp_get_num_procs();
	const int MaxThreads = nCPUCore;
	cout << "                             nCPUCore[" << nCPUCore << "],MaxThreads[" << MaxThreads << "]" << endl;
	int omp_threads = nthread > MaxThreads ? MaxThreads : nthread;
	if (omp_threads < 1) {
		omp_threads = MaxThreads;
	}
	int nScreenThread = omp_threads;
	if (nScreenThread > NBUS) {
		nScreenThread = NBUS;
	}
//#pragma omp parallel for num_threads(nScreenThread)	//并行计算for循环，暂时封上
	for (int i = sNum; i < eNum; i++)
	{
		size_t index = 0;
		size_t inbbs = 0, iwbbs = 0;
		double *ZZ1 = NULL;
		inbbs = bus[i].ibs;//内部计算母线号
		if (NBUS == 1) {
			ZZ1 = ZZ1_bs;
		}
		else {
			double *NodeRX = GB1.RetIMp();
			ZZ1 = NodeRX + nrow * 2 * inbbs;
		}
		//故障点短路电流计算
		double vbs = bus[i].vc, ai = bus[i].ac, vbase = bus[i].vbase;
		double IBase = wbase*_GH3 / vbase * 1000;//故障点短路电流
		MYCOMPLEX ZZ1_ii;//自阻抗
		MYCOMPLEX zf;
		index = inbbs * 2;
		ZZ1_ii.r = ZZ1[index];
		ZZ1_ii.i = ZZ1[index + 1];
		zf.r = rf;
		zf.i = xf;
		scc3_out[i].v_fp.r = vbs*cos(ai);
		scc3_out[i].v_fp.i = vbs*sin(ai);
		scc3_out[i].if3 = scc3_out[i].v_fp / (ZZ1_ii + zf);
		//scc3_out[i].clear();	//测试封上，否则会清空结果
		scc3_out[i].i_fp = scc3_out[i].if3.r*scc3_out[i].if3.r + scc3_out[i].if3.i*scc3_out[i].if3.i;
		scc3_out[i].i_fp = sqrt(scc3_out[i].i_fp);
		scc3_out[i].mva_fp = scc3_out[i].i_fp*bus[i].vc*wbase;
		scc3_out[i].i_fp = scc3_out[i].i_fp*IBase;
		scc3_out[i].VBus_scc[i].ibs = bus[i].id;
		scc3_out[i].VBus_scc[i].name = bus[i].name;
		scc3_out[i].VBus_scc[i].kvvl = bus[i].vbase;
		scc3_out[i].VBus_scc[i].iisland = iisland;
		scc3_out[i].VBus_scc[i].bsindex = i;
		scc3_out[i].VBus_scc[i].v = scc3_out[i].v_fp;
		//其他任意节点电压计算
		MYCOMPLEX ZZ1_ij;
		//scc3_out[i].VBus_scc.clear();	//测试封上，否则会清空结果
		//scc3_out[i].VBus_scc.resize(nbus);	//测试封上，否则会清空结果
		for (int ii = 0; ii < nrow; ii++)//节点导纳编号
		{
			iwbbs = inb_iwb[ii];
			if (ii == inbbs) {
				continue;
			}
			index = ii * 2;
			ZZ1_ij.r = ZZ1[index];
			ZZ1_ij.i = ZZ1[index + 1];
			scc3_out[i].VBus_scc[iwbbs].name = bus[iwbbs].name;
			scc3_out[i].VBus_scc[iwbbs].ibs = bus[iwbbs].id;
			scc3_out[i].VBus_scc[iwbbs].iisland = iisland;
			scc3_out[i].VBus_scc[iwbbs].bsindex = iwbbs;
			scc3_out[i].VBus_scc[iwbbs].kvvl = bus[iwbbs].vbase;
			scc3_out[i].VBus_scc[iwbbs].v_pre.r = bus[iwbbs].vc*cos(bus[iwbbs].ac);
			scc3_out[i].VBus_scc[iwbbs].v_pre.i = bus[iwbbs].vc*sin(bus[iwbbs].ac);
			scc3_out[i].VBus_scc[iwbbs].v = scc3_out[i].VBus_scc[iwbbs].v_pre - ZZ1_ij*scc3_out[i].if3;
			scc3_out[i].VBus_scc[iwbbs].vthree = scc3_out[i].VBus_scc[iwbbs].v.r*scc3_out[i].VBus_scc[iwbbs].v.r + scc3_out[i].VBus_scc[iwbbs].v.i + scc3_out[i].VBus_scc[iwbbs].v.i;
			scc3_out[i].VBus_scc[iwbbs].vthree = sqrt(scc3_out[i].VBus_scc[iwbbs].vthree);
		}
		//任意支路短路电流
		MYCOMPLEX zij;
		MYCOMPLEX Iij_bch;
		size_t iwb = 0, jwb = 0;
		scc3_out[i].IBranch_scc.clear();
		scc3_out[i].IBranch_scc.resize(nbranch);
		for (int ii = 0; ii < nbranch; ii++)
		{
			zij.r = branch[ii].r;
			zij.i = branch[ii].x;
			iwb = branch[ii].i;
			jwb = branch[ii].j;
			scc3_out[i].IBranch_scc[ii].ibranch = ii;
			if (branch[ii].branchtype == BRANCHTYPE_LN) {//线路
				Iij_bch.r = 0;
				Iij_bch.i = branch[ii].bch2;
			}
			else if (branch[ii].branchtype == BRANCHTYPE_LNS) {//串补
				Iij_bch.r = 0;
				Iij_bch.i = 0;
			}
			else if (branch[ii].branchtype == BRANCHTYPE_XF) {//变压器
				Iij_bch.r = 0;
				Iij_bch.i = 0;
			}
			scc3_out[i].IBranch_scc[ii].i3_i = (scc3_out[i].VBus_scc[iwb].v - scc3_out[i].VBus_scc[jwb].v) / zij + scc3_out[i].VBus_scc[iwb].v*Iij_bch;
			scc3_out[i].IBranch_scc[ii].i3_j = (scc3_out[i].VBus_scc[jwb].v - scc3_out[i].VBus_scc[iwb].v) / zij + scc3_out[i].VBus_scc[jwb].v*Iij_bch;
			scc3_out[i].IBranch_scc[ii].i_i = sqrt(scc3_out[i].IBranch_scc[ii].i3_i.r*scc3_out[i].IBranch_scc[ii].i3_i.r + scc3_out[i].IBranch_scc[ii].i3_i.i*scc3_out[i].IBranch_scc[ii].i3_i.i);
			scc3_out[i].IBranch_scc[ii].i_j = sqrt(scc3_out[i].IBranch_scc[ii].i3_j.r*scc3_out[i].IBranch_scc[ii].i3_j.r + scc3_out[i].IBranch_scc[ii].i3_j.i*scc3_out[i].IBranch_scc[ii].i3_j.i);
			vbase = bus[iwb].vbase;
			IBase = wbase*_GH3 / vbase * 1000;
			scc3_out[i].IBranch_scc[ii].i_i *= IBase;
			scc3_out[i].IBranch_scc[ii].i_j *= IBase;
		}
	}

	return 1;
}
int CIeeeSCCBase::SinglePhaseGroundShortCircuit(int bs_index, int nthread, eABCPhase sphase, double rf, double xf)//单相接地短路
{
	if (nbus == 1) {
		ErrorMessage << "就1个节点，开玩笑？" << endl;
		return -1;
	}
	int NBUS = 0;//扫描母线总数
	int sNum = 0;//扫描起始位置
	int eNum = 0;//扫描结束位置
	double *ZZ1_bs = NULL;
	double *ZZ2_bs = NULL;
	double *ZZ0_bs = NULL;
	const size_t nrow = GB1.GetIDim();//各序节点阻抗矩阵行数
	if (bs_index < 0) {//扫描所有母线单相接地短路故障
		NBUS = nbus;
		sNum = 0;
		for (int i = 0; i < NBUS; i++)	//初始化节点计算结果vector和支路计算结果vector的维数
		{
			scc1_out[i].VBus_scc.resize(NBUS);
			scc1_out[i].IBranch_scc.resize(nbranch);
		}
		if (GB1.Inversion() < 0) {
			ErrorMessage << "矩阵求逆失败！" << endl;
			return -1;
		}
		if (GB2.Inversion() < 0) {
			ErrorMessage << "矩阵求逆失败！" << endl;
			return -1;
		}
		if (GB0.Inversion() < 0) {
			ErrorMessage << "矩阵求逆失败！" << endl;
			return -1;
		}
	}
	else if (bs_index>nbus) {
		ErrorMessage << "指定母线下表[" << bs_index << "]无效！" << endl;
		return -1;
	}
	else {
		NBUS = 1;
		sNum = bs_index;
		scc1_out[sNum].VBus_scc.resize(nrow);	//初始化节点计算结果vector的维数
		scc1_out[sNum].IBranch_scc.resize(nbranch);	//初始化支路计算结果vector的维数
		if (GB1.Inversion(bus[bs_index].ibs, &ZZ1_bs) < 0) {
			ErrorMessage << "求解bus[" << bs_index << "][bs|" << bus[bs_index].id << "](" << bus[bs_index].name << ")正序自阻抗互阻抗失败!" << endl;
			return -1;
		}
		if (GB2.Inversion(bus[bs_index].ibs, &ZZ2_bs) < 0) {
			ErrorMessage << "求解bus[" << bs_index << "][bs|" << bus[bs_index].id << "](" << bus[bs_index].name << ")负序自阻抗互阻抗失败!" << endl;
			return -1;
		}
		if (GB0.Inversion(bus[bs_index].ibs, &ZZ0_bs) < 0) {
			ErrorMessage << "求解bus[" << bs_index << "][bs|" << bus[bs_index].id << "](" << bus[bs_index].name << ")零序自阻抗互阻抗失败!" << endl;
			return -1;
		}
	}
	eNum = sNum + NBUS;
	int *inb_iwb = GetNBI_WBI();
	//设置并行线程
	const int nCPUCore = omp_get_num_procs();
	const int MaxThreads = nCPUCore;
	cout << "                             nCPUCore[" << nCPUCore << "],MaxThreads[" << MaxThreads << "]" << endl;
	int omp_threads = nthread > MaxThreads ? MaxThreads : nthread;
	if (omp_threads < 1) {
		omp_threads = MaxThreads;
	}
	int nScreenThread = omp_threads;
	if (nScreenThread > NBUS) {
		nScreenThread = NBUS;
	}
#pragma omp parallel for num_threads(nScreenThread)
	for (int i = sNum; i < eNum; i++)
	{
		size_t index = 0;
		size_t inbbs = 0, iwbbs = 0;
		double *ZZ1 = NULL;
		double *ZZ2 = NULL;
		double *ZZ0 = NULL;
		inbbs = bus[i].ibs;//内部计算母线号
		if (NBUS == 1) {
			ZZ1 = ZZ1_bs;
			ZZ2 = ZZ2_bs;
			ZZ0 = ZZ0_bs;
		}
		else {
			double *NodeRX = GB1.RetIMp();
			ZZ1 = NodeRX + nrow * 2 * inbbs; 
			NodeRX = GB2.RetIMp();
			ZZ2 = NodeRX + nrow * 2 * inbbs; 
			NodeRX = GB0.RetIMp();
			ZZ0 = NodeRX + nrow * 2 * inbbs;
		}
		//故障点短路电流计算
		double vbs = bus[i].vc, ai = bus[i].ac, vbase = bus[i].vbase;
		double IBase = wbase*_GH3 / vbase * 1000;//故障点短路电流
		MYCOMPLEX ZZ1_ii, ZZ2_ii, ZZ0_ii;//各序自阻抗
		MYCOMPLEX zf;//短路阻抗
		index = inbbs * 2;
		ZZ1_ii.r = ZZ1[index];
		ZZ1_ii.i = ZZ1[index + 1];
		ZZ2_ii.r = ZZ2[index];
		ZZ2_ii.i = ZZ2[index + 1];
		ZZ0_ii.r = ZZ0[index];
		ZZ0_ii.i = ZZ0[index + 1];
		zf.r = rf;
		zf.i = xf;
		UnsymmetricalShortcircuit_OUT *scc1 = &scc1_out[i];
		scc1->v_fp.r = vbs*cos(ai);
		scc1->v_fp.i = vbs*sin(ai);
		scc1->i1_fp = scc1->v_fp / (ZZ1_ii + ZZ2_ii + ZZ0_ii + zf);
		scc1->i2_fp = scc1->i1_fp;
		scc1->i0_fp = scc1->i1_fp;
		scc1->v1_fp = scc1->v_fp - scc1->i1_fp*ZZ1_ii;
		scc1->v2_fp = -1 * scc1->i2_fp*ZZ2_ii;
		scc1->v0_fp = -1 * scc1->i0_fp*ZZ0_ii;

		//其他任意节点电压计算
		MYCOMPLEX ZZ1_ij;
		MYCOMPLEX ZZ2_ij;
		MYCOMPLEX ZZ0_ij;
		scc1->VBus_scc.clear();
		scc1->VBus_scc.resize(nbus);
		for (int ii = 0; ii < nrow; ii++)//节点导纳编号
		{
			iwbbs = inb_iwb[ii];
			if (ii == inbbs) {
				continue;
			}
			index = ii * 2;
			ZZ1_ij.r = ZZ1[index];
			ZZ1_ij.i = ZZ1[index + 1];
			ZZ2_ij.r = ZZ2[index];
			ZZ2_ij.i = ZZ2[index + 1];
			ZZ0_ij.r = ZZ0[index];
			ZZ0_ij.i = ZZ0[index + 1];
			scc1->VBus_scc[iwbbs].name = bus[iwbbs].name;
			scc1->VBus_scc[iwbbs].ibs = bus[iwbbs].id;
			scc1->VBus_scc[iwbbs].iisland = iisland;
			scc1->VBus_scc[iwbbs].bsindex = iwbbs;
			scc1->VBus_scc[iwbbs].kvvl = bus[iwbbs].vbase;
			scc1->VBus_scc[iwbbs].v_pre.r = bus[iwbbs].vc*cos(bus[iwbbs].ac);
			scc1->VBus_scc[iwbbs].v_pre.i = bus[iwbbs].vc*sin(bus[iwbbs].ac);
			scc1->VBus_scc[iwbbs].v1 = scc1->VBus_scc[iwbbs].v_pre - ZZ1_ij*scc1->i1_fp;
			scc1->VBus_scc[iwbbs].v2 = -1 * ZZ2_ij*scc1->i2_fp;
			scc1->VBus_scc[iwbbs].v0 = -1 * ZZ0_ij*scc1->i0_fp;
		}
		//任意支路短路电流
		MYCOMPLEX zij1, zij2, zij0;
		MYCOMPLEX Iij1_bch, Iij2_bch, Iij0_bch;
		size_t iwb = 0, jwb = 0;
		scc1->IBranch_scc.clear();
		scc1->IBranch_scc.resize(nbranch);
		for (int ii = 0; ii < nbranch; ii++)
		{
			zij1.r = branch[ii].r;
			zij1.i = branch[ii].x;
			zij2.r = branch[ii].r;
			zij2.i = branch[ii].x;
			zij0.r = branch[ii].r0;
			zij0.i = branch[ii].x0;
			iwb = branch[ii].i;
			jwb = branch[ii].j;
			scc1->IBranch_scc[ii].ibranch = ii;
			if (branch[ii].branchtype == BRANCHTYPE_LN) {//线路
				Iij1_bch.r = 0;
				Iij1_bch.i = branch[ii].bch2;
				Iij2_bch.r = 0;
				Iij2_bch.i = branch[ii].bch2;
				Iij0_bch.r = 0;
				Iij0_bch.i = branch[ii].bch0_2;
			}
			else if (branch[ii].branchtype == BRANCHTYPE_LNS) {//串补
				Iij1_bch.r = 0;
				Iij1_bch.i = 0;
				Iij2_bch.r = 0;
				Iij2_bch.i = 0;
				Iij0_bch.r = 0;
				Iij0_bch.i = 0;
			}
			else if (branch[ii].branchtype == BRANCHTYPE_XF) {//变压器
				Iij1_bch.r = 0;
				Iij1_bch.i = 0;
				Iij2_bch.r = 0;
				Iij2_bch.i = 0;
				Iij0_bch.r = 0;
				Iij0_bch.i = 0;
			}
			scc1->IBranch_scc[ii].i1_i = (scc1->VBus_scc[iwb].v1 - scc1->VBus_scc[jwb].v1) / zij1 + scc1->VBus_scc[iwb].v1*Iij1_bch;
			scc1->IBranch_scc[ii].i2_i = (scc1->VBus_scc[iwb].v2 - scc1->VBus_scc[jwb].v2) / zij2 + scc1->VBus_scc[iwb].v2*Iij2_bch;
			scc1->IBranch_scc[ii].i0_i = (scc1->VBus_scc[iwb].v0 - scc1->VBus_scc[jwb].v0) / zij0 + scc1->VBus_scc[iwb].v0*Iij0_bch;
			scc1->IBranch_scc[ii].i1_j = (scc1->VBus_scc[jwb].v1 - scc1->VBus_scc[iwb].v1) / zij1 + scc1->VBus_scc[jwb].v1*Iij1_bch;
			scc1->IBranch_scc[ii].i2_j = (scc1->VBus_scc[jwb].v2 - scc1->VBus_scc[iwb].v2) / zij2 + scc1->VBus_scc[jwb].v2*Iij2_bch;
			scc1->IBranch_scc[ii].i0_j = (scc1->VBus_scc[jwb].v0 - scc1->VBus_scc[iwb].v0) / zij0 + scc1->VBus_scc[jwb].v0*Iij0_bch;
			//scc1->IBranch_scc[ii].i_i = sqrt(scc1->IBranch_scc[ii].i3_i.r*scc1->IBranch_scc[ii].i3_i.r + scc1->IBranch_scc[ii].i3_i.i*scc1->IBranch_scc[ii].i3_i.i);
			//scc1->IBranch_scc[ii].i_j = sqrt(scc1->IBranch_scc[ii].i3_j.r*scc1->IBranch_scc[ii].i3_j.r + scc1->IBranch_scc[ii].i3_j.i*scc1->IBranch_scc[ii].i3_j.i);
			//vbase = bus[iwb].vbase;
			//IBase = wbase*_GH3 / vbase * 1000;
			//scc1->IBranch_scc[ii].i_i *= IBase;
			//scc1->IBranch_scc[ii].i_j *= IBase;
		}
	}

	return 1;
}
int CIeeeSCCBase::TwoPhaseShortCircuit(int bs_index, eABCPhase sphase, double rf, double xf)//两相相间短路
{
	//不用算零序，零序都为0
	double *ZZ1 = NULL;
	double *ZZ2 = NULL;
	double *ZZ0 = NULL;
	if (GB1.Inversion(bus[bs_index].ibs, &ZZ1) < 0) {
		ErrorMessage << "求解bus[" << bs_index << "][bs|" << bus[bs_index].id << "](" << bus[bs_index].name << ")正序自阻抗互阻抗失败!" << endl;
		return -1;
	}
	if (GB2.Inversion(bus[bs_index].ibs, &ZZ2) < 0) {
		ErrorMessage << "求解bus[" << bs_index << "][bs|" << bus[bs_index].id << "](" << bus[bs_index].name << ")负序自阻抗互阻抗失败!" << endl;
		return -1;
	}
	//if (GB0.Inversion(bus[bs_index].ibs, &ZZ0) < 0) {
	//	ErrorMessage << "求解bus[" << bs_index << "][bs|" << bus[bs_index].id << "](" << bus[bs_index].name << ")零序自阻抗互阻抗失败!" << endl;
	//	return -1;
	//}
	//故障点短路电流计算
	double vbs = bus[bs_index].vc, ai = bus[bs_index].ac, vbase = bus[bs_index].vbase;
	double IBase = wbase*_GH3 / vbase * 1000;//故障点短路电流
	MYCOMPLEX ZZ1_ii, ZZ2_ii;// , ZZ0_ii;//各序自阻抗
	MYCOMPLEX zf;//短路阻抗
	size_t index = 0;
	size_t inbbs = 0, iwbbs = 0;
	inbbs = bus[bs_index].ibs;//内部计算母线号
	index = inbbs * 2;
	ZZ1_ii.r = ZZ1[index];
	ZZ1_ii.i = ZZ1[index + 1];
	ZZ2_ii.r = ZZ2[index];
	ZZ2_ii.i = ZZ2[index + 1];
	//ZZ0_ii.r = ZZ0[index];
	//ZZ0_ii.i = ZZ0[index + 1];
	zf.r = rf;
	zf.i = xf;
	scc2.clear();
	scc2.v_fp.r = vbs*cos(ai);
	scc2.v_fp.i = vbs*sin(ai);
	scc2.i1_fp = scc2.v_fp / (ZZ1_ii + ZZ2_ii + zf);
	scc2.i2_fp = -1 * scc2.i1_fp;
	scc2.i0_fp.r = 0;
	scc2.i0_fp.i = 0;
	scc2.v1_fp = scc2.v_fp - scc2.i1_fp*ZZ1_ii;
	scc2.v2_fp = -1 * scc2.i2_fp*ZZ2_ii;
	scc2.v0_fp.r = 0;
	scc2.v0_fp.i = 0;

	const size_t nrow = GB1.GetIDim();//各序节点阻抗矩阵行数
	int *inb_iwb = GetNBI_WBI();
	//其他任意节点电压计算
	MYCOMPLEX ZZ1_ij;
	MYCOMPLEX ZZ2_ij;
	//MYCOMPLEX ZZ0_ij;
	scc2.VBus_scc.clear();
	scc2.VBus_scc.resize(nbus);
	for (int ii = 0; ii < nrow; ii++)//节点导纳编号
	{
		iwbbs = inb_iwb[ii];
		if (ii == inbbs) {
			continue;
		}
		index = ii * 2;
		ZZ1_ij.r = ZZ1[index];
		ZZ1_ij.i = ZZ1[index + 1];
		ZZ2_ij.r = ZZ2[index];
		ZZ2_ij.i = ZZ2[index + 1];
		//ZZ0_ij.r = ZZ0[index];
		//ZZ0_ij.i = ZZ0[index + 1];
		scc2.VBus_scc[iwbbs].name = bus[iwbbs].name;
		scc2.VBus_scc[iwbbs].ibs = bus[iwbbs].id;
		scc2.VBus_scc[iwbbs].iisland = iisland;
		scc2.VBus_scc[iwbbs].bsindex = iwbbs;
		scc2.VBus_scc[iwbbs].kvvl = bus[iwbbs].vbase;
		scc2.VBus_scc[iwbbs].v_pre.r = bus[iwbbs].vc*cos(bus[iwbbs].ac);
		scc2.VBus_scc[iwbbs].v_pre.i = bus[iwbbs].vc*sin(bus[iwbbs].ac);
		scc2.VBus_scc[iwbbs].v1 = scc2.VBus_scc[iwbbs].v_pre - ZZ1_ij*scc2.i1_fp;
		scc2.VBus_scc[iwbbs].v2 = -1 * ZZ2_ij*scc2.i2_fp;
		//scc2.VBus_scc[iwbbs].v0 = -1 * ZZ0_ij*scc2.i0_fp;//=0
		scc2.VBus_scc[iwbbs].v0.r = 0;
		scc2.VBus_scc[iwbbs].v0.i = 0;
	}
	//任意支路短路电流
	MYCOMPLEX zij1, zij2;// , zij0;
	MYCOMPLEX Iij1_bch, Iij2_bch;// , Iij0_bch;
	size_t iwb = 0, jwb = 0;
	scc2.IBranch_scc.clear();
	scc2.IBranch_scc.resize(nbranch);
	for (int ii = 0; ii < nbranch; ii++)
	{
		zij1.r = branch[ii].r;
		zij1.i = branch[ii].x;
		zij2.r = branch[ii].r;
		zij2.i = branch[ii].x;
		//zij0.r = branch[ii].r0;
		//zij0.i = branch[ii].x0;
		iwb = branch[ii].i;
		jwb = branch[ii].j;
		scc2.IBranch_scc[ii].ibranch = ii;
		if (branch[ii].branchtype == BRANCHTYPE_LN) {//线路
			Iij1_bch.r = 0;
			Iij1_bch.i = branch[ii].bch2;
			Iij2_bch.r = 0;
			Iij2_bch.i = branch[ii].bch2;
			//Iij0_bch.r = 0;
			//Iij0_bch.i = branch[ii].bch0_2;
		}
		else if (branch[ii].branchtype == BRANCHTYPE_LNS) {//串补
			Iij1_bch.r = 0;
			Iij1_bch.i = 0;
			Iij2_bch.r = 0;
			Iij2_bch.i = 0;
			//Iij0_bch.r = 0;
			//Iij0_bch.i = 0;
		}
		else if (branch[ii].branchtype == BRANCHTYPE_XF) {//变压器
			Iij1_bch.r = 0;
			Iij1_bch.i = 0;
			Iij2_bch.r = 0;
			Iij2_bch.i = 0;
			//Iij0_bch.r = 0;
			//Iij0_bch.i = 0;
		}
		scc2.IBranch_scc[ii].i1_i = (scc2.VBus_scc[iwb].v1 - scc2.VBus_scc[jwb].v1) / zij1 + scc2.VBus_scc[iwb].v1*Iij1_bch;
		scc2.IBranch_scc[ii].i2_i = (scc2.VBus_scc[iwb].v2 - scc2.VBus_scc[jwb].v2) / zij2 + scc2.VBus_scc[iwb].v2*Iij2_bch;
		//scc2.IBranch_scc[ii].i0_i = (scc2.VBus_scc[iwb].v0 - scc2.VBus_scc[jwb].v0) / zij0 + scc2.VBus_scc[iwb].v0*Iij0_bch;//=0
		scc2.IBranch_scc[ii].i0_i.r = 0;
		scc2.IBranch_scc[ii].i0_i.i = 0;
		scc2.IBranch_scc[ii].i1_j = (scc2.VBus_scc[jwb].v1 - scc2.VBus_scc[iwb].v1) / zij1 + scc2.VBus_scc[jwb].v1*Iij1_bch;
		scc2.IBranch_scc[ii].i2_j = (scc2.VBus_scc[jwb].v2 - scc2.VBus_scc[iwb].v2) / zij2 + scc2.VBus_scc[jwb].v2*Iij2_bch;
		//scc2.IBranch_scc[ii].i0_j = (scc2.VBus_scc[jwb].v0 - scc2.VBus_scc[iwb].v0) / zij0 + scc2.VBus_scc[jwb].v0*Iij0_bch;//=0
		scc2.IBranch_scc[ii].i0_j.r = 0;
		scc2.IBranch_scc[ii].i0_j.i = 0;
	}

	return 1;
}
int CIeeeSCCBase::TwoPhaseGroundShortCircuit(int bs_index, eABCPhase sphase, double rf, double xf)//两相接地短路
{
	double *ZZ1 = NULL;
	double *ZZ2 = NULL;
	double *ZZ0 = NULL;
	if (GB1.Inversion(bus[bs_index].ibs, &ZZ1) < 0) {
		ErrorMessage << "求解bus[" << bs_index << "][bs|" << bus[bs_index].id << "](" << bus[bs_index].name << ")正序自阻抗互阻抗失败!" << endl;
		return -1;
	}
	if (GB2.Inversion(bus[bs_index].ibs, &ZZ2) < 0) {
		ErrorMessage << "求解bus[" << bs_index << "][bs|" << bus[bs_index].id << "](" << bus[bs_index].name << ")负序自阻抗互阻抗失败!" << endl;
		return -1;
	}
	if (GB0.Inversion(bus[bs_index].ibs, &ZZ0) < 0) {
		ErrorMessage << "求解bus[" << bs_index << "][bs|" << bus[bs_index].id << "](" << bus[bs_index].name << ")零序自阻抗互阻抗失败!" << endl;
		return -1;
	}
	//故障点短路电流计算
	double vbs = bus[bs_index].vc, ai = bus[bs_index].ac, vbase = bus[bs_index].vbase;
	double IBase = wbase*_GH3 / vbase * 1000;//故障点短路电流
	MYCOMPLEX ZZ1_ii, ZZ2_ii, ZZ0_ii;//各序自阻抗
	MYCOMPLEX zf;//短路阻抗
	size_t index = 0;
	size_t inbbs = 0, iwbbs = 0;
	inbbs = bus[bs_index].ibs;//内部计算母线号
	index = inbbs * 2;
	ZZ1_ii.r = ZZ1[index];
	ZZ1_ii.i = ZZ1[index + 1];
	ZZ2_ii.r = ZZ2[index];
	ZZ2_ii.i = ZZ2[index + 1];
	ZZ0_ii.r = ZZ0[index];
	ZZ0_ii.i = ZZ0[index + 1];
	zf.r = rf;
	zf.i = xf;
	scc2.clear();
	scc2.v_fp.r = vbs*cos(ai);
	scc2.v_fp.i = vbs*sin(ai);
	scc2.i1_fp = scc2.v_fp / (ZZ1_ii + ((ZZ2_ii*(ZZ0_ii + 3 * zf)) / (ZZ2_ii + ZZ0_ii + 3 * zf)));
	scc2.i2_fp = -1 * scc2.i1_fp*(ZZ0_ii + 3 * zf)/(ZZ2_ii + ZZ0_ii + 3 * zf);
	scc2.i0_fp = -1 * scc2.i1_fp*ZZ2_ii / (ZZ2_ii + ZZ0_ii + 3 * zf);
	scc2.v1_fp = scc2.v_fp - scc2.i1_fp*ZZ1_ii;
	scc2.v2_fp = -1 * scc2.i2_fp*ZZ2_ii;
	scc2.v0_fp = -1 * scc2.i0_fp*ZZ0_ii;

	const size_t nrow = GB1.GetIDim();//各序节点阻抗矩阵行数
	int *inb_iwb = GetNBI_WBI();
	//其他任意节点电压计算
	MYCOMPLEX ZZ1_ij;
	MYCOMPLEX ZZ2_ij;
	MYCOMPLEX ZZ0_ij;
	scc2.VBus_scc.clear();
	scc2.VBus_scc.resize(nbus);
	for (int ii = 0; ii < nrow; ii++)//节点导纳编号
	{
		iwbbs = inb_iwb[ii];
		if (ii == inbbs) {
			continue;
		}
		index = ii * 2;
		ZZ1_ij.r = ZZ1[index];
		ZZ1_ij.i = ZZ1[index + 1];
		ZZ2_ij.r = ZZ2[index];
		ZZ2_ij.i = ZZ2[index + 1];
		ZZ0_ij.r = ZZ0[index];
		ZZ0_ij.i = ZZ0[index + 1];
		scc2.VBus_scc[iwbbs].name = bus[iwbbs].name;
		scc2.VBus_scc[iwbbs].ibs = bus[iwbbs].id;
		scc2.VBus_scc[iwbbs].iisland = iisland;
		scc2.VBus_scc[iwbbs].bsindex = iwbbs;
		scc2.VBus_scc[iwbbs].kvvl = bus[iwbbs].vbase;
		scc2.VBus_scc[iwbbs].v_pre.r = bus[iwbbs].vc*cos(bus[iwbbs].ac);
		scc2.VBus_scc[iwbbs].v_pre.i = bus[iwbbs].vc*sin(bus[iwbbs].ac);
		scc2.VBus_scc[iwbbs].v1 = scc2.VBus_scc[iwbbs].v_pre - ZZ1_ij*scc2.i1_fp;
		scc2.VBus_scc[iwbbs].v2 = -1 * ZZ2_ij*scc2.i2_fp;
		scc2.VBus_scc[iwbbs].v0 = -1 * ZZ0_ij*scc2.i0_fp;
	}
	//任意支路短路电流
	MYCOMPLEX zij1, zij2, zij0;
	MYCOMPLEX Iij1_bch, Iij2_bch, Iij0_bch;
	size_t iwb = 0, jwb = 0;
	scc2.IBranch_scc.clear();
	scc2.IBranch_scc.resize(nbranch);
	for (int ii = 0; ii < nbranch; ii++)
	{
		zij1.r = branch[ii].r;
		zij1.i = branch[ii].x;
		zij2.r = branch[ii].r;
		zij2.i = branch[ii].x;
		zij0.r = branch[ii].r0;
		zij0.i = branch[ii].x0;
		iwb = branch[ii].i;
		jwb = branch[ii].j;
		scc2.IBranch_scc[ii].ibranch = ii;
		if (branch[ii].branchtype == BRANCHTYPE_LN) {//线路
			Iij1_bch.r = 0;
			Iij1_bch.i = branch[ii].bch2;
			Iij2_bch.r = 0;
			Iij2_bch.i = branch[ii].bch2;
			Iij0_bch.r = 0;
			Iij0_bch.i = branch[ii].bch0_2;
		}
		else if (branch[ii].branchtype == BRANCHTYPE_LNS) {//串补
			Iij1_bch.r = 0;
			Iij1_bch.i = 0;
			Iij2_bch.r = 0;
			Iij2_bch.i = 0;
			Iij0_bch.r = 0;
			Iij0_bch.i = 0;
		}
		else if (branch[ii].branchtype == BRANCHTYPE_XF) {//变压器
			Iij1_bch.r = 0;
			Iij1_bch.i = 0;
			Iij2_bch.r = 0;
			Iij2_bch.i = 0;
			Iij0_bch.r = 0;
			Iij0_bch.i = 0;
		}
		scc2.IBranch_scc[ii].i1_i = (scc2.VBus_scc[iwb].v1 - scc2.VBus_scc[jwb].v1) / zij1 + scc2.VBus_scc[iwb].v1*Iij1_bch;
		scc2.IBranch_scc[ii].i2_i = (scc2.VBus_scc[iwb].v2 - scc2.VBus_scc[jwb].v2) / zij2 + scc2.VBus_scc[iwb].v2*Iij2_bch;
		scc2.IBranch_scc[ii].i0_i = (scc2.VBus_scc[iwb].v0 - scc2.VBus_scc[jwb].v0) / zij0 + scc2.VBus_scc[iwb].v0*Iij0_bch;//=0
		scc2.IBranch_scc[ii].i1_j = (scc2.VBus_scc[jwb].v1 - scc2.VBus_scc[iwb].v1) / zij1 + scc2.VBus_scc[jwb].v1*Iij1_bch;
		scc2.IBranch_scc[ii].i2_j = (scc2.VBus_scc[jwb].v2 - scc2.VBus_scc[iwb].v2) / zij2 + scc2.VBus_scc[jwb].v2*Iij2_bch;
		scc2.IBranch_scc[ii].i0_j = (scc2.VBus_scc[jwb].v0 - scc2.VBus_scc[iwb].v0) / zij0 + scc2.VBus_scc[jwb].v0*Iij0_bch;//=0
	}

	return 1;
}
int CIeeeSCCBase::SinglePhaseBreak(int ibs_index, int jbs_index, eABCPhase sphase)//一相断线
{

	if (nbus == 1) {
		ErrorMessage << "就1个节点，开玩笑？" << endl;
		return -1;
	}
	int NBUS = 0;//扫描母线总数
	int sNum = 0;//扫描起始位置
	int eNum = 0;//扫描结束位置
	double *ZZ1_bs = NULL;
	double *ZZ2_bs = NULL;
	double *ZZ0_bs = NULL;
	const size_t nrow = GB1.GetIDim();//各序节点阻抗矩阵行数
	if (ibs_index > nbus) {
		ErrorMessage << "指定母线下表[" << ibs_index << "]无效！" << endl;
		return -1;
	}
	else {
		NBUS = 1;
		sNum = ibs_index;
		m_scc1Break_out.VBus_scc.resize(nrow);	//初始化节点计算结果vector的维数
		m_scc1Break_out.IBranch_scc.resize(nbranch);	//初始化支路计算结果vector的维数
		if (GB1.Inversion(bus[ibs_index].ibs, &ZZ1_bs) < 0) {
			ErrorMessage << "求解bus[" << ibs_index << "][bs|" << bus[ibs_index].id << "](" << bus[ibs_index].name << ")正序自阻抗互阻抗失败!" << endl;
			return -1;
		}
		if (GB2.Inversion(bus[ibs_index].ibs, &ZZ2_bs) < 0) {
			ErrorMessage << "求解bus[" << ibs_index << "][bs|" << bus[ibs_index].id << "](" << bus[ibs_index].name << ")负序自阻抗互阻抗失败!" << endl;
			return -1;
		}
		if (GB0.Inversion(bus[ibs_index].ibs, &ZZ0_bs) < 0) {
			ErrorMessage << "求解bus[" << ibs_index << "][bs|" << bus[ibs_index].id << "](" << bus[ibs_index].name << ")零序自阻抗互阻抗失败!" << endl;
			return -1;
		}
	}
	eNum = sNum + NBUS;
	int *inb_iwb = GetNBI_WBI();
	
	for (int i = sNum; i < eNum; i++)
	{
		size_t index = 0;
		size_t inbbs = 0, iwbbs = 0;
		double *ZZ1 = NULL;
		double *ZZ2 = NULL;
		double *ZZ0 = NULL;
		inbbs = bus[i].ibs;//内部计算母线号
		if (NBUS == 1) {
			ZZ1 = ZZ1_bs;
			ZZ2 = ZZ2_bs;
			ZZ0 = ZZ0_bs;
		}
		//故障点短路电流计算
		double vbs = bus[i].vc, ai = bus[i].ac, vbase = bus[i].vbase;
		double IBase = wbase * _GH3 / vbase * 1000;//故障点短路电流
		MYCOMPLEX ZZ1_ii, ZZ2_ii, ZZ0_ii;//各序自阻抗
		MYCOMPLEX zf;//短路阻抗
		index = inbbs * 2;
		ZZ1_ii.r = ZZ1[index];
		ZZ1_ii.i = ZZ1[index + 1];
		ZZ2_ii.r = ZZ2[index];
		ZZ2_ii.i = ZZ2[index + 1];
		ZZ0_ii.r = ZZ0[index];
		ZZ0_ii.i = ZZ0[index + 1];
		UnsymmetricalShortcircuit_OUT *scc1 = &m_scc1Break_out;
		scc1->v_fp.r = vbs * cos(ai);
		scc1->v_fp.i = vbs * sin(ai);
		scc1->i1_fp = scc1->v_fp / (ZZ1_ii + ZZ2_ii + ZZ0_ii + zf);
		scc1->i2_fp = scc1->i1_fp;
		scc1->i0_fp = scc1->i1_fp;
		scc1->v1_fp = scc1->v_fp - scc1->i1_fp*ZZ1_ii;
		scc1->v2_fp = -1 * scc1->i2_fp*ZZ2_ii;
		scc1->v0_fp = -1 * scc1->i0_fp*ZZ0_ii;

		//其他任意节点电压计算
		MYCOMPLEX ZZ1_ij;
		MYCOMPLEX ZZ2_ij;
		MYCOMPLEX ZZ0_ij;
		scc1->VBus_scc.clear();
		scc1->VBus_scc.resize(nbus);
		for (int ii = 0; ii < nrow; ii++)//节点导纳编号
		{
			iwbbs = inb_iwb[ii];
			if (ii == inbbs) {
				continue;
			}
			index = ii * 2;
			ZZ1_ij.r = ZZ1[index];
			ZZ1_ij.i = ZZ1[index + 1];
			ZZ2_ij.r = ZZ2[index];
			ZZ2_ij.i = ZZ2[index + 1];
			ZZ0_ij.r = ZZ0[index];
			ZZ0_ij.i = ZZ0[index + 1];
			scc1->VBus_scc[iwbbs].name = bus[iwbbs].name;
			scc1->VBus_scc[iwbbs].ibs = bus[iwbbs].id;
			scc1->VBus_scc[iwbbs].iisland = iisland;
			scc1->VBus_scc[iwbbs].bsindex = iwbbs;
			scc1->VBus_scc[iwbbs].kvvl = bus[iwbbs].vbase;
			scc1->VBus_scc[iwbbs].v_pre.r = bus[iwbbs].vc*cos(bus[iwbbs].ac);
			scc1->VBus_scc[iwbbs].v_pre.i = bus[iwbbs].vc*sin(bus[iwbbs].ac);
			scc1->VBus_scc[iwbbs].v1 = scc1->VBus_scc[iwbbs].v_pre - ZZ1_ij * scc1->i1_fp;
			scc1->VBus_scc[iwbbs].v2 = -1 * ZZ2_ij*scc1->i2_fp;
			scc1->VBus_scc[iwbbs].v0 = -1 * ZZ0_ij*scc1->i0_fp;
		}
		//任意支路短路电流
		MYCOMPLEX zij1, zij2, zij0;
		MYCOMPLEX Iij1_bch, Iij2_bch, Iij0_bch;
		size_t iwb = 0, jwb = 0;
		scc1->IBranch_scc.clear();
		scc1->IBranch_scc.resize(nbranch);
		for (int ii = 0; ii < nbranch; ii++)
		{
			zij1.r = branch[ii].r;
			zij1.i = branch[ii].x;
			zij2.r = branch[ii].r;
			zij2.i = branch[ii].x;
			zij0.r = branch[ii].r0;
			zij0.i = branch[ii].x0;
			iwb = branch[ii].i;
			jwb = branch[ii].j;
			scc1->IBranch_scc[ii].ibranch = ii;
			if (branch[ii].branchtype == BRANCHTYPE_LN) {//线路
				Iij1_bch.r = 0;
				Iij1_bch.i = branch[ii].bch2;
				Iij2_bch.r = 0;
				Iij2_bch.i = branch[ii].bch2;
				Iij0_bch.r = 0;
				Iij0_bch.i = branch[ii].bch0_2;
			}
			else if (branch[ii].branchtype == BRANCHTYPE_LNS) {//串补
				Iij1_bch.r = 0;
				Iij1_bch.i = 0;
				Iij2_bch.r = 0;
				Iij2_bch.i = 0;
				Iij0_bch.r = 0;
				Iij0_bch.i = 0;
			}
			else if (branch[ii].branchtype == BRANCHTYPE_XF) {//变压器
				Iij1_bch.r = 0;
				Iij1_bch.i = 0;
				Iij2_bch.r = 0;
				Iij2_bch.i = 0;
				Iij0_bch.r = 0;
				Iij0_bch.i = 0;
			}
			scc1->IBranch_scc[ii].i1_i = (scc1->VBus_scc[iwb].v1 - scc1->VBus_scc[jwb].v1) / zij1 + scc1->VBus_scc[iwb].v1*Iij1_bch;
			scc1->IBranch_scc[ii].i2_i = (scc1->VBus_scc[iwb].v2 - scc1->VBus_scc[jwb].v2) / zij2 + scc1->VBus_scc[iwb].v2*Iij2_bch;
			scc1->IBranch_scc[ii].i0_i = (scc1->VBus_scc[iwb].v0 - scc1->VBus_scc[jwb].v0) / zij0 + scc1->VBus_scc[iwb].v0*Iij0_bch;
			scc1->IBranch_scc[ii].i1_j = (scc1->VBus_scc[jwb].v1 - scc1->VBus_scc[iwb].v1) / zij1 + scc1->VBus_scc[jwb].v1*Iij1_bch;
			scc1->IBranch_scc[ii].i2_j = (scc1->VBus_scc[jwb].v2 - scc1->VBus_scc[iwb].v2) / zij2 + scc1->VBus_scc[jwb].v2*Iij2_bch;
			scc1->IBranch_scc[ii].i0_j = (scc1->VBus_scc[jwb].v0 - scc1->VBus_scc[iwb].v0) / zij0 + scc1->VBus_scc[jwb].v0*Iij0_bch;
			//scc1->IBranch_scc[ii].i_i = sqrt(scc1->IBranch_scc[ii].i3_i.r*scc1->IBranch_scc[ii].i3_i.r + scc1->IBranch_scc[ii].i3_i.i*scc1->IBranch_scc[ii].i3_i.i);
			//scc1->IBranch_scc[ii].i_j = sqrt(scc1->IBranch_scc[ii].i3_j.r*scc1->IBranch_scc[ii].i3_j.r + scc1->IBranch_scc[ii].i3_j.i*scc1->IBranch_scc[ii].i3_j.i);
			//vbase = bus[iwb].vbase;
			//IBase = wbase*_GH3 / vbase * 1000;
			//scc1->IBranch_scc[ii].i_i *= IBase;
			//scc1->IBranch_scc[ii].i_j *= IBase;
		}
	}

	return 1;
}
int CIeeeSCCBase::TwoPhaseBreak(int ibs_index, int jbs_index, eABCPhase sphase)//两相断线
{
	return 1;
}
int CIeeeSCCBase::SCCCalc(eUnsymmetricalShortcircuitType scc_type, int ibs_index, int jbs_index)
{
	sSCCType(scc_type);
	if (scc_type == eThreePhaseShortcircuit) {
		scc3_out.clear();
		scc3_out.resize(nbus);
		GetSCCGB(1);
		ThreePhaseShortcircuit(1);
	}
	else if (scc_type == eSinglePhaseGroundShortCircuit) {
		scc1_out.clear();
		scc1_out.resize(nbus);
		GetSCCGB(0);	//计算零负序导纳矩阵
		SinglePhaseGroundShortCircuit(1);
	}
	else if (scc_type == eTwoPhaseShortCircuit) {
		scc1_out.clear();
		scc1_out.resize(nbus);
		GetSCCGB(0);	//计算零负序导纳矩阵
		TwoPhaseShortCircuit(1, eAPhase);
	}
	else if (scc_type == eTwoPhaseGroundShortCircuit) {
		scc1_out.clear();
		scc1_out.resize(nbus);
		GetSCCGB(0);	//计算零负序导纳矩阵
		TwoPhaseGroundShortCircuit(1, eAPhase);
	}
	printGB(1, 1, 1);

	return 1;
}

CIeeeSCC::CIeeeSCC()
{
	nbus = 0;
	nbranch = 0;
	NodeRX = NULL;
	If3 = NULL;
	Id3 = NULL;
}
CIeeeSCC::~CIeeeSCC()
{
	if (If3)delete[]If3;
	if (Id3 != NULL && nbus > 0) {
		for (int i = 0; i < nbus; i++)
		{
			if (Id3[i])delete Id3[i];
		}
		delete[]Id3;
	}
}
int CIeeeSCC::GetNodeS3GB()
{
	if (!IsReady()) {
		AlarmMessage << "Ieee模型未准备完毕！" << endl;
		return -1;
	}
#ifdef _DEBUG
	TIMING jsks("CIeeeSCC::GetNodeS3GB");
#endif
	int it = 0;
	double g = 0, b = 0, bch2 = 0, t = 0;
	int i = 0, j = 0;
	LISTELE * iigb = NULL;
	LISTELE * jjgb = NULL;
	LISTELE * ijgb = NULL;
	AlarmMessage << "节点导纳矩阵维数(" << iNAllACNode << "*" << iNAllACNode << ")" << endl;
	S3GB.init(iNAllACNode);

	for (it = 0; it < nbranch; it++)
	{
		if (branch[it].branchtype == BRANCHTYPE_LN) {//线路
			i = branch[it].i;
			j = branch[it].j;
			i = bus[i].ibs;//转换内部节点
			j = bus[j].ibs;//转换内部节点
			g = branch[it].g;
			b = branch[it].b;
			bch2 = branch[it].bch2;///10000;//1/2线路电纳
									  //branch[it].bch2 = 0;
			iigb = S3GB.FindOrAdd(i, i);//S3GB为对称矩阵，仅保存上三角
			jjgb = S3GB.FindOrAdd(j, j);
			ijgb = S3GB.FindOrAdd(i, j);
			//节点i,j自阻抗
			//节点i,j自阻抗
			iigb->real += g;
			iigb->imag += (b + bch2);
			jjgb->real += g;
			jjgb->imag += (b + bch2);
			//i,j互阻抗
			ijgb->real -= g;
			ijgb->imag -= b;
#ifdef _DEBUG
			if (S3GB.EleIsNanOrInf(iigb) || S3GB.EleIsNanOrInf(jjgb) || S3GB.EleIsNanOrInf(ijgb)) {
				AlarmMessage << "ln_branch[" << it << "],name(" << branch[it].name << "),g(" << branch[it].g << "),b(" << branch[it].b << ")" << endl;
				return -1;
			}
#endif
		}
		else if (branch[it].branchtype == BRANCHTYPE_LNS) {//串补
			i = branch[it].i;
			j = branch[it].j;
			i = bus[i].ibs;//转换内部节点
			j = bus[j].ibs;//转换内部节点
			b = branch[it].b;
			iigb = S3GB.FindOrAdd(i, i);//S3GB为对称矩阵，仅保存上三角
			jjgb = S3GB.FindOrAdd(j, j);
			ijgb = S3GB.FindOrAdd(i, j);
			//节点i,j自阻抗
			//节点i,j自阻抗
			iigb->imag += b;
			jjgb->imag += b;
			//i,j互阻抗
			ijgb->imag -= b;
#ifdef _DEBUG
			if (S3GB.EleIsNanOrInf(iigb) || S3GB.EleIsNanOrInf(jjgb) || S3GB.EleIsNanOrInf(ijgb)) {
				AlarmMessage << "lns_branch[" << it << "],name(" << branch[it].name << "),b(" << branch[it].b << ")" << endl;
				return -1;
			}
#endif
		}
		else if (branch[it].branchtype == BRANCHTYPE_XF) {//变压器支路//t在一端，i端
			g = branch[it].g;
			b = branch[it].b;
			t = branch[it].t;
			if (t<0.6 || t>1.4) {
				AlarmMessage << "xf_branch[" << it << "],name(" << branch[it].name << "),g(" << g << "),b(" << b << "),t(" << t << "),变比异常！" << endl;
				t = 1;
			}
			i = branch[it].i;//高压侧，变压器参数折算到高压侧，因此，i侧为标准侧，j侧为非标准差（状态估计书中的叫法）?????????
			j = branch[it].j;//低压侧
			i = bus[i].ibs;//转换内部节点
			j = bus[j].ibs;//转换内部节点
			iigb = S3GB.FindOrAdd(i, i);
			jjgb = S3GB.FindOrAdd(j, j);
			ijgb = S3GB.FindOrAdd(i, j);
			//i端为非标准侧，j端标准侧//节点i,j自阻抗(pwr_new中也这么处理的，fd_yy函数中tj=1)
			//branch[it].t = 1;
			iigb->real += (g / (t*t));
			iigb->imag += (b / (t*t));
			jjgb->real += (g);
			jjgb->imag += (b);
			//i,j互阻抗
			ijgb->real -= (g / t);
			ijgb->imag -= (b / t);
#ifdef _DEBUG
			if (S3GB.EleIsNanOrInf(iigb) || S3GB.EleIsNanOrInf(jjgb) || S3GB.EleIsNanOrInf(ijgb)) {
				AlarmMessage << "xf_branch[" << it << "],name(" << branch[it].name << "),g(" << branch[it].g << "),b(" << branch[it].b << "),t(" << branch[it].t << ")" << endl;
				return -1;
			}
#endif
		}
		else {
			AlarmMessage << "支路：branch[" << it << "](" << branch[it].name << ")->branchtype=" << branch[it].branchtype << "未处理！！！" << endl;
		}
	}
	//cout<<endl<<endl;

	for (it = 0; it < nbus; it++)//并联电容器、电抗器、负荷、机组
	{
		//if (DCBus_DCIsland[it] > 0)continue;//直流母线
		if (bus[it].nodetype == NODE_DCBS)continue;//直流母线，正常是不会有并联容抗器的
		if (bus[it].nacbranch == 0)continue;//无交流支路，即该节点为换流变侧交流节点，通过换流器与主网有连接，但换流变开关断开
		i = bus[it].ibs;
		iigb = S3GB.FindOrAdd(i, i);
		//并联电容器/电抗器
		iigb->real += bus[it].blg;
		iigb->imag += bus[it].blb;
		if (bus[it].nld > 0 || bus_ndccnv[it] > 0) {//负荷、换流器采用恒阻抗模型
			bus[it].ldg = ((bus[it].ldp)*_wbase) / (bus[it].vc*bus[it].vc);
			bus[it].ldb = (-(bus[it].ldq)*_wbase) / (bus[it].vc*bus[it].vc);
			//bus[it].ldg = ((bus[it].ldp - bus[it].dccnvp)*_wbase) / (bus[it].vc*bus[it].vc);
			//bus[it].ldb = (-(bus[it].ldq - bus[it].dccnvq)*_wbase) / (bus[it].vc*bus[it].vc);
			//AlarmMessage << "[bs|" << bus[it].id << "](" << bus[it].ibs << ")(" << bus[it].ldg << "," << bus[it].ldb << ")(" << bus[it].ldp << "," << bus[it].ldq << ")" << endl;
			//if (fabs(bus[it].dccnvp) > 0.0) {
			//	AlarmMessage << "[bs|" << bus[it].id << "](" << it << "|" << bus[it].ibs << ")dccnv(" << bus[it].dccnvp << "," << bus[it].dccnvq << ")ld(" << bus[it].ldp << "," << bus[it].ldq << ")换流变侧母线！" << endl;
			//	AlarmMessage << "(" << bus[it].ldg << "," << bus[it].ldb << ")" << "(" << (bus[it].ldp*_wbase) / (bus[it].vc*bus[it].vc) << "," << (-bus[it].ldq*_wbase) / (bus[it].vc*bus[it].vc) << ")" << endl;
			//}
			iigb->real += bus[it].ldg;
			iigb->imag += bus[it].ldb;
		}
		if (bus[it].nun > 0) {//机组次暂态电抗，没参数则采用恒阻抗模型			
			if (fabs(bus[it].ungd) > EPSINON) {//机组
				//AlarmMessage << "[bs|" << bus[it].id << "].ungd=" << bus[it].ungd << "!" << endl;
				iigb->imag += bus[it].ungd;
			}
			else {
				//if (fabs(bus[it].unp) > 1.0 || fabs(bus[it].unq) > 1.0)
				{
					//AlarmMessage << "匹配不上的机组按负负荷处理！" << endl;
					iigb->real += (-bus[it].unp*_wbase) / (bus[it].vc*bus[it].vc);
					iigb->imag += (bus[it].unq*_wbase) / (bus[it].vc*bus[it].vc);
				}
			}
		}

#ifdef _DEBUG
		if (S3GB.EleIsNanOrInf(iigb)) {
			return -1;
		}
#endif
	}
	S3GB.SetIsReady(true);
	return 1;
}
int CIeeeSCC::SCCModelInit(IeeeToolBase *MSrc)
{
	if (IeeeToolBaseCopy(MSrc) < 0) {
		AlarmMessage << "拷贝模型IeeeToolBaseCopy错误!" << endl;
		return -1;
	}
	if (GetNodeS3GB() < 0) {
		AlarmMessage << "获取系统母线三项短路时节点导纳矩阵（负荷采用恒定导纳模型）失败！" << endl;
		return -1;
	}
	if (S3GB_KLU.InitBySCM_Column_KLU(&S3GB, eComplexNum2, 0) < 0) {
		AlarmMessage << "KLU节点导纳矩阵初始化错误！！！" << endl;
		return -1;
	}
	if (S3GB_KLU.Inversion() < 0) {
		AlarmMessage << "导纳矩阵求逆错误！！！" << endl;
		return -1;
	}
	NodeRX = S3GB_KLU.RetIMp();
	if (NodeRX == NULL) {
		AlarmMessage << "获取系统短路时节点阻抗矩阵（负荷采用恒定阻抗模型）失败，节点导纳矩阵求逆失败！" << endl;
		return -1;
	}
	If3 = new double[nbus];
	Id3 = new double*[nbus];
	for (int i = 0; i < nbus; i++)
	{
		Id3[i] = NULL;
	}
	return 1;
}
int CIeeeSCC::ScanBSSCC_3(int ibs)
{
	TIMING ksjs("ScanBSSCC_3");
	int ret = 0;
	int NBUS = 0;//扫描母线总数
	int sNum = 0;//扫描起始位置
	int eNum = 0;//扫描结束位置
	if (ibs < 0) {//扫描所有母线三项短路
		NBUS = nbus;
		sNum = 0;
	}
	else {
		NBUS = 1;
		sNum = ibs;
	}
	eNum = sNum + NBUS;

	int *inb_iwb = GetNBI_WBI();
	size_t nrow = S3GB.RIDim();
#pragma omp parallel for
	for (int i = sNum; i < eNum; i++)
	{
		size_t index = 0;
		size_t inb = 0, jnb = 0, inbbs = 0, iwbbs = 0;
		double Rii = 0.0, Xii = 0.0;//节点阻抗矩阵自阻抗
		double Rij = 0.0, Xij = 0.0;//节点阻抗矩阵互阻抗
		double Zii = 0.0, Zij = 0.0;//自(互)阻抗的模
		double vbs = 0.0, vbase = 0.0;
		double ai = 0;
		MYCOMPLEX *Vi = NULL;
		MYCOMPLEX VVector;//电压向量
		MYCOMPLEX IF3, RX, _Id3;// , RXij, RXji;
		double IBase = 0.0;//故障点短路电流

		Id3[i] = new double[nbranch];//支路短路电流
		inbbs = bus[i].ibs;//内部计算母线号
		vbs = bus[i].vc;
		ai = bus[i].ac;
		VVector.r = -vbs*cos(ai);
		VVector.i = -vbs*sin(ai);
		vbase = bus[i].vbase;
		IBase = wbase*_GH3 / vbase * 1000;
		index = inbbs + inbbs*nrow;
		Rii = CREAL(NodeRX, index);
		Xii = CIMAG(NodeRX, index);
		RX.r = Rii;
		RX.i = Xii;
		Zii = Rii*Rii + Xii*Xii;
		Zii = sqrt(Zii);
		IF3 = VVector / RX;//复数之比的模等于模之比
		//vbs = 1;
		If3[i] = -vbs / Zii;
		//If3 = -_GH3*vbs / Zii;
		//If3 = If3*IBase;
		//cout <<"If="<< If3[i] << " (节点电压=" << vbs << "|v=" << bus[i].v << ") 短路容量=" << If3[i]* vbs * 100 <<"("<< 100 * 1.05*1.05 / Zii <<") 短路电流="<< If3[i]* IBase << endl;// " " << sqrt(IF3.r*IF3.r + IF3.i*IF3.i) << endl;// "(" << GH3*If3 / IBase << ")" << " 短路容量：" << wbase * (vbs*vbs / Zii) << " " << If3*vbase*vbs*GH3*0.001 << endl;
		/*cout << "If=" << If3[i] << " (节点电压=" << vbs << "|v=" << bus[i].v << ") 短路容量=" << If3[i] * vbs * 100 << " 短路电流=" << If3[i] * IBase << " IBase=" << IBase << endl;// " " << sqrt(IF3.r*IF3.r + IF3.i*IF3.i) << endl;// "(" << GH3*If3 / IBase << ")" << " 短路容量：" << wbase * (vbs*vbs / Zii) << " " << If3*vbase*vbs*GH3*0.001 << endl;
		cout << "r=" << RX.r << " x=" << RX.i << endl;*/
		
		Vi = new MYCOMPLEX[nrow];
		for (size_t j = 0; j < nrow; j++)
		{
			index = j + inbbs*nrow;
			RX.r = greal(NodeRX, index);
			RX.i = gimag(NodeRX, index);
			Vi[j] = RX*IF3;
			iwbbs = inb_iwb[j];
			vbs = bus[iwbbs].vc;
			ai = bus[iwbbs].ac;
			VVector.r = vbs*cos(ai);
			VVector.i = vbs*sin(ai);
			Vi[j] = VVector - Vi[j];//故障分量+故障前电压
			//cout << Vi[j] << endl;
		}
		for (size_t j = 0; j < nbranch; j++)
		{
			RX.r = branch[j].r;
			RX.i = branch[j].x;
			inb = branch[j].i;
			inb = bus[inb].ibs;
			jnb = branch[j].j;
			jnb = bus[jnb].ibs;
			_Id3 = (Vi[inb] - Vi[jnb]) / RX;
			Id3[i][j] = sqrt(_Id3.r*_Id3.r + _Id3.i*_Id3.i);
			vbase = bus[inb].vbase;
			IBase = wbase*_GH3 / vbase * 1000;
			//cout << "支路(" << branch[j].name << ")短路电流标幺值(" << Id3[i][j] << ")" << _Id3 << endl;
		}

		delete[]Vi;
	}

	return 1;
}

int ShortCircuitRatioCalc::ReadSetMiescr(string fname)
{
	ifstream flog;
	string line;
	SetMiescr SMTemp;
	SCRClear();
	flog.open(fname.c_str());
	if (!flog.is_open())
	{
		cout << "打开文件" << fname << "失败！" << endl;
		exit(0);
	}
	while (!flog.eof()) {
		line.clear();
		getline(flog, line, '\n');
		if (line.find_first_of('\r') != string::npos) {
			line = line.substr(0, line.find_first_of('\r'));
		}
		if (line.erase(0, line.find_first_not_of("	 ")).empty())continue;
		line >> SMTemp.miescrname >> SMTemp.stname >> SMTemp.kv >> SMTemp.mva;
		vSetMiescr.push_back(SMTemp);
	}
	flog.close();
	lv_miescr = vSetMiescr.size();
	nSCRBus = lv_miescr;
	DCCNVName = new string[nSCRBus];
	fscc = new double[nSCRBus];
	fscc2 = new double[nSCRBus];
	MIESCR2 = new double[nSCRBus];
	MIESCR = new double[nSCRBus];
	QF = new double[nSCRBus];
	MIIF = new double*[nSCRBus];
	IF3 = new double[nSCRBus];
	iRX = new double*[nSCRBus];
	RX = new MYCOMPLEX*[nSCRBus];
	iSCRBus = new int[nSCRBus];
	iNBBS = new int[nSCRBus];
	IsDCCNV = new bool[nSCRBus];

	for (int j = 0; j < nSCRBus; j++)
	{
		fscc[j] = 0.0;
		fscc2[j] = 0.0;
		MIESCR2[j] = 0.0;
		MIESCR[j] = 0.0;
		QF[j] = 0.0;
		MIIF[j] = NULL;
		IF3[j] = 0.0;
		iRX[j] = NULL;
		RX[j] = NULL;
		iSCRBus[j] = -1;
		iNBBS[j] = -1;
		IsDCCNV[j] = false;
	}
	for (int j = 0; j < nSCRBus; j++)
	{
		for (int i = 0; i < scc.nbus; i++)
		{
			//if ((strstr(scc.bus[i].name, vSetMiescr[j].stname.c_str()) != NULL))
			//{
			//	AlarmMessage << scc.bus[i].name << " " << vSetMiescr[j].stname << " " << scc.bus[i].vbase << endl;
			//}
			if ((strstr(scc.bus[i].name, vSetMiescr[j].stname.c_str()) != NULL) && (fabs(scc.bus[i].vbase - vSetMiescr[j].kv) < EPSINON)) {
				iSCRBus[j] = i;
				iNBBS[j] = scc.bus[i].ibs;
				DCCNVName[j] = vSetMiescr[j].miescrname;
				AlarmMessage << scc.bus[i].name << " " << DCCNVName[j] << " iSCRBus[" << j << "]=" << iSCRBus[j] << " iNBBS[" << j << "]=" << iNBBS[j] << endl;
				break;
			}
		}
	}

	return 1;
}
void ShortCircuitRatioCalc::SCRClear()
{
	lv_miescr = 0;
	vSetMiescr.clear();
	if (DCCNVName)delete[]DCCNVName;//直流名称
	DCCNVName = NULL;
	if (iNBBS)delete[]iNBBS;//内部节点号集合
	iNBBS = NULL;
	if (iSCRBus)delete[]iSCRBus;//外部节点号集合，即bus下标
	iSCRBus = NULL;
	if (IsDCCNV)delete[]IsDCCNV;//节点换流器是否投运
	IsDCCNV = NULL;
	if (IF3)delete[]IF3;//节点短路电流
	IF3 = NULL;
	if (iRX) {
		for (size_t i = 0; i < nSCRBus; i++)
		{
			if (iRX[i]) {
				delete[]iRX[i];
			}
		}
		delete[]iRX;
		iRX = NULL;
	}
	if (RX) {
		for (size_t i = 0; i < nSCRBus; i++)
		{
			if (RX[i]) {
				delete[]RX[i];
			}
		}
		delete[]RX;
		RX = NULL;
	}
	if (MIIF) {
		for (size_t i = 0; i < nSCRBus; i++)
		{
			if (MIIF[i]) {
				delete[]MIIF[i];
			}
		}
		delete[]MIIF;
		MIIF = NULL;
	}
	if (QF)delete[]QF;//无功补偿
	QF = NULL;
	if (MIESCR)delete[]MIESCR;//短路比100*v/z
	MIESCR = NULL;
	if (MIESCR2)delete[]MIESCR2;//短路比v*v/z
	MIESCR2 = NULL;
	if (fscc)delete[]fscc;//短路容量100*v/z
	fscc = NULL;
	if (fscc2)delete[]fscc2;//短路容量v*v/z
	fscc2 = NULL;
	nSCRBus = 0;
	if (dccnvp)delete[]dccnvp;
	dccnvp = NULL;
	if (dccnvq)delete[]dccnvq;
	dccnvq = NULL;
}
int ShortCircuitRatioCalc::SCRInit(IeeeToolBase *MSrc)
{
	if (scc.IeeeToolBaseCopy(MSrc) < 0) {
		AlarmMessage << "拷贝模型IeeeToolBaseCopy错误!" << endl;
		return -1;
	}
	if (scc.GetNodeS3GB() < 0) {
		AlarmMessage << "生成零阻抗节点三项短路节点导纳矩阵GetNodeS3GB错误!" << endl;
		return -1;
	}
	if (scc.S3GB_KLU.InitBySCM_Column_KLU(&scc.S3GB, eComplexNum2, 0) < 0) {
		AlarmMessage << "KLU节点导纳矩阵初始化错误！！！" << endl;
		return -1;
	}
	return 1;
}
int ShortCircuitRatioCalc::ShortCircuitRatio()
{
	if (nSCRBus <= 0) {
		AlarmMessage << "输入计算母线数为" << nSCRBus << ",短路比计算ShortCircuitRatio失败！" << endl;
		return -1;
	}
	string stemp;
	//nSCRBus = nbs;
	size_t iSCRBranch = 0;
	dccnvp = new double[nSCRBus];
	dccnvq = new double[nSCRBus];
	AlarmMessage << endl << "+++++++++++++++++++++++++++++++++计算短路比的节点支路如下+++++++++++++++++++++++++++++++++" << endl;
	for (size_t i = 0; i < nSCRBus; i++)
	{
		if (iSCRBus[i] < 0 || iSCRBus[i] >= scc.nbus) {
			AlarmMessage << "输入计算母线逻辑号错误" << iSCRBus[i] << "[0," << scc.nbus << "],短路比计算ShortCircuitRatio失败！" << endl;
			return -1;
		}
		AlarmMessage << "输入计算母线逻辑号=" << iSCRBus[i] << endl;
		fscc[i] = 0;
		fscc2[i] = 0;
		MIESCR2[i] = 0;
		MIESCR[i] = 0;
		QF[i] = 0;
		MIIF[i] = NULL;
		IF3[i] = 0;
		iRX[i] = NULL;
		RX[i] = NULL;
		dccnvp[i] = 0.0;
		dccnvq[i] = 0.0;
		for (size_t j = 0; j < scc.bus_ndccnv[iSCRBus[i]]; j++)
		{
			dccnvp[i] += scc.dccnv[scc.bus_dccnv[iSCRBus[i]][j]].pac;
			dccnvq[i] += scc.dccnv[scc.bus_dccnv[iSCRBus[i]][j]].qac;
		}
		//if (fabs(bus[iSCRBus[i]].dccnvp) < 10.0) {// && fabs(bus[iSCRBus[i]].dccnvq) < 10.0) {
		if (fabs(dccnvp[i]) < 10.0) {// && fabs(bus[iSCRBus[i]].dccnvq) < 10.0) {
			stemp = "换流器停运";
			IsDCCNV[i] = false;
		}
		else {
			stemp = "换流器投运";
			IsDCCNV[i] = true;
		}
		AlarmMessage << "bs(" << iSCRBus[i] << ")[bs|" << scc.bus[iSCRBus[i]].id << "]" << "(nb|" << scc.bus[iSCRBus[i]].ibs << ")(v=" << scc.bus[iSCRBus[i]].vc << "," << scc.bus[iSCRBus[i]].v << ")该节点[" << stemp << "],功率[" << dccnvp[i] << "," << dccnvq[i] << "]" << endl;
		for (size_t j = 0; j < scc.bus_nbranch[iSCRBus[i]]; j++)
		{
			iSCRBranch = scc.bus_branch[iSCRBus[i]][j];
			AlarmMessage << "		branch[" << iSCRBranch << "](" << scc.branch[iSCRBranch].name << ")" << endl;
		}
	}
	AlarmMessage << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

	SuiteSparse_long index = 0;
	size_t index2 = 0;
	double Zij = 0;
	MYCOMPLEX cMIIF;

	//去掉换流站的所有无功补偿
	for (size_t i = 0; i < nSCRBus; i++)
	{
		if (!IsDCCNV[i])continue;
		index = -1;
		for (SuiteSparse_long ll = scc.S3GB_KLU.Ap[iNBBS[i]]; ll < scc.S3GB_KLU.Ap[iNBBS[i] + 1]; ll++)//循环列数据
		{
			if (scc.S3GB_KLU.Ai[ll] == iNBBS[i]) {//行号等于列号，则为对角元
				index = 2 * ll;
				break;
			}
		}
		if (index == -1) {
			AlarmMessage << "不可能吧，节点导纳对角元为0？" << endl;
			return -1;
		}
		//去掉该节点无功补偿装置(不包括线路高抗)
		//AlarmMessage << "去掉该节点[bs|"<< scc.bus[iSCRBus[i]].id <<"](" << scc.bus[iSCRBus[i]].name << ")无功补偿装置(不包括线路高抗)!" << endl;
		scc.S3GB_KLU.Ax[index + 1] -= scc.bus[iSCRBus[i]].bsb;
		//AlarmMessage << scc.S3GB_KLU.Ax[index + 1] << endl;
	}

	for (size_t i = 0; i < nSCRBus; i++)
	{
		if (!IsDCCNV[i])continue;
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//index = -1;
		//for (SuiteSparse_long ll = S3GB_KLU.Ap[iNBBS[i]]; ll < S3GB_KLU.Ap[iNBBS[i] + 1]; ll++)//循环列数据
		//{
		//	if (S3GB_KLU.Ai[ll] == iNBBS[i]) {//行号等于列号，则为对角元
		//		index = 2 * ll;
		//		break;
		//	}
		//}
		//if (index == -1) {
		//	AlarmMessage << "不可能吧，节点导纳对角元为0？" << endl;
		//	return -1;
		//}
		//去掉该节点无功补偿装置
		//S3GB_KLU.Ax[index + 1] -= bus[iSCRBus[i]].bsb;
		//S3GB_KLU.Ax[index] -= bus[iSCRBus[i]].blg;
		//S3GB_KLU.Ax[index + 1] -= bus[iSCRBus[i]].blb;
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		if (scc.S3GB_KLU.Inversion(iNBBS[i], &iRX[i]) < 0) {
			AlarmMessage << "求解(" << iSCRBus[i] << ")[bs|" << scc.bus[iSCRBus[i]].id << "]自阻抗互阻抗失败!" << endl;
			return -1;
		}
		RX[i] = new MYCOMPLEX[nSCRBus];
		MIIF[i] = new double[nSCRBus];
		for (size_t ii = 0; ii < nSCRBus; ii++)
		{
			if (!IsDCCNV[ii]) {
				RX[i][ii].r = 0;
				RX[i][ii].i = 0;
			}
			else {
				index2 = iNBBS[ii] * 2;
				RX[i][ii].r = iRX[i][index2];
				RX[i][ii].i = iRX[i][index2 + 1];
			}
		}
		Zij = 0.0;
		Zij = sqrt(RX[i][i].r*RX[i][i].r + RX[i][i].i*RX[i][i].i);
		if (Zij < EPSINON) {
			AlarmMessage << "(" << iSCRBus[i] << ")[bs|" << scc.bus[iSCRBus[i]].id << "]自阻抗模为(" << Zij << ")，不可能吧！" << endl;
			return -1;
		}
		QF[i] = scc.wbase *scc.bus[iSCRBus[i]].blb * scc.bus[iSCRBus[i]].vc*scc.bus[iSCRBus[i]].vc;
		IF3[i] = scc.bus[iSCRBus[i]].vc / Zij;
		//AlarmMessage << "短路电流标幺值：" << IF3[i] << endl;
		for (size_t ii = 0; ii < nSCRBus; ii++)
		{
			if (IsDCCNV[ii]) {
				cMIIF = RX[i][ii] / RX[i][i];
				MIIF[i][ii] = sqrt(cMIIF.r*cMIIF.r + cMIIF.i*cMIIF.i);
			}
			else {
				MIIF[i][ii] = 0;
			}
		}

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//恢复该节点无功补偿装置
		//S3GB_KLU.Ax[index + 1] += bus[iSCRBus[i]].bsb;
		//S3GB_KLU.Ax[index] += bus[iSCRBus[i]].blg;
		//S3GB_KLU.Ax[index + 1] += bus[iSCRBus[i]].blb;
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	}

	for (int i = 0; i < nSCRBus; i++)
	{
		fscc[i] = 0;
		fscc2[i] = 0;
		if (IsDCCNV[i]) {
			fscc[i] = IF3[i] * 100;
			fscc2[i] = IF3[i] * scc.bus[iSCRBus[i]].vc * 100;
		}
	}
	double dtemp = 0;
	for (size_t i = 0; i < nSCRBus; i++)
	{
		if (IsDCCNV[i]) {
			MIESCR[i] = 0;
			MIESCR2[i] = 0;
			dtemp = 0;
			for (size_t ii = 0; ii < nSCRBus; ii++)
			{
				if (IsDCCNV[ii]) {
					if (i != ii) {
						dtemp += (fabs(dccnvp[ii])*MIIF[i][ii]);
					}
				}
			}
			MIESCR[i] = (fscc[i] - fabs(QF[i])) / (fabs(dccnvp[i]) + dtemp);
			MIESCR2[i] = (fscc2[i] - fabs(QF[i])) / (fabs(dccnvp[i]) + dtemp);
		}
	}

	return 1;
}
int ShortCircuitRatioCalc::PrintMIESCR()
{
	//结果输出
	ofstream log;
	log.open("短路比结果.csv");
	if (!log.is_open()) {
		AlarmMessage << "打开文件（短路比结果.csv）失败!" << endl;
		return 1;
	}
	log << "节点导纳矩阵（节点阻抗矩阵）维数（" << scc.S3GB.idim << "*" << scc.S3GB.jdim << "）" << endl;
	log << "<阻抗表>" << endl;
	log << "直流名称,";
	for (size_t i = 0; i < nSCRBus; i++)
	{
		if (IsDCCNV[i]) {
			if (i < (nSCRBus - 1))log << DCCNVName[i] << ",";
			else log << DCCNVName[i] << endl;
		}
		else {
			if (i == (nSCRBus - 1)) {
				log << endl;
			}
		}
	}
	for (size_t i = 0; i < nSCRBus; i++)
	{
		if (IsDCCNV[i]) {
			log << DCCNVName[i] << ",";
			for (size_t ii = 0; ii < nSCRBus; ii++)
			{
				if (IsDCCNV[ii]) {
					if (ii < (nSCRBus - 1))
						log << RX[i][ii].r << "|" << RX[i][ii].i << ",";
					else
						log << RX[i][ii].r << "|" << RX[i][ii].i << endl;
				}
				else {
					if (ii == (nSCRBus - 1)) {
						log << endl;
					}
				}
			}
		}
	}
	log << endl << endl;

	log << "<短路比结果表>" << endl;
	log << "表,电压,电压量测,电阻,电抗,短路容量(短路电流标幺值*100),短路容量(v*v/z),pdc,qf,短路比(短路电流标幺值*100),短路比(v*v/z)" << endl;
	for (size_t i = 0; i < nSCRBus; i++)
	{
		if (IsDCCNV[i]) {
			log << DCCNVName[i] << "," << scc.bus[iSCRBus[i]].vc << "," << scc.bus[iSCRBus[i]].v << "," << RX[i][i].r << "," << RX[i][i].i << "," << fscc[i] << "," << fscc2[i] << "," << dccnvp[i] << "," << QF[i] << "," << MIESCR[i] << "," << MIESCR2[i] << endl;
		}
	}
	log << endl << endl;

	log << "<影响影子>" << endl;
	log << "直流名称,";
	for (size_t i = 0; i < nSCRBus; i++)
	{
		if (IsDCCNV[i]) {
			if (i < (nSCRBus - 1))log << DCCNVName[i] << ",";
			else log << DCCNVName[i] << endl;
		}
		else {
			if (i == (nSCRBus - 1)) {
				log << endl;
			}
		}
	}
	for (size_t i = 0; i < nSCRBus; i++)
	{
		if (IsDCCNV[i]) {
			log << DCCNVName[i] << ",";
			for (size_t ii = 0; ii < nSCRBus; ii++)
			{
				if (IsDCCNV[ii]) {
					if (ii < (nSCRBus - 1))
						log << MIIF[i][ii] << ",";
					else
						log << MIIF[i][ii] << endl;
				}
				else {
					if (ii == (nSCRBus - 1)) {
						log << endl;
					}
				}
			}
		}
	}
	log << endl << endl;
	log.close();
	return 1;
}
int CIeeeSCCBase::printGB(int gb1_flag, int gb2_flag, int gb0_flag)
{
	string fileName = "C:\\Users\\sunbo\\Desktop\\1.PSASource\\debug\\shortCircuitType";
	if (gb1_flag == 1)
	{
		//fileName += "_gb1.txt";
		int line = gb1.RIDim();
		int row = gb1.RJDim();
		gb1.printM(fileName + "_gb1.txt");
		GB1.PrintSM(fileName + "_GB1_klu.txt");
	}
	if (gb2_flag == 1)
	{
		//fileName += "_gb2.txt";
		int line = gb2.RIDim();
		int row = gb2.RJDim();
		gb2.printM(fileName + "_gb2.txt");
		GB2.PrintSM(fileName + "_GB2_klu.txt");
	}
	if (gb0_flag == 1)
	{
		//fileName += "_gb0.txt";
		int line = gb0.RIDim();
		int row = gb0.RJDim();
		gb0.printM(fileName + "_gb0.txt");
		GB0.PrintSM(fileName + "_GB0_klu.txt");
	}
	return 0;
}