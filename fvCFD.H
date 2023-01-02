#ifndef fvCFD_H
#define fvCFD_H

#include "parRun.H"       //并行计算

#include "Time.H"         //时间项
#include "fvMesh.H"       //有限体积法离散所需要的网格信息
#include "fvc.H"          //方程求解显式项，向量形式
#include "fvMatrices.H"   //方程求解隐式项，矩阵形式
#include "fvm.H"          //与上一行相关，返回矩阵
#include "linear.H"       //线性，与插值相关（如gauss linear）
#include "uniformDimensionedFields.H"   //
#include "calculatedFvPatchFields.H"    //边界条件场，类型为calculated，一般边界设置随时间变化时采用该格式
#include "extrapolatedCalculatedFvPatchFields.H"    //边界条件类型，外插
#include "fixedValueFvPatchFields.H"    //边界条件，固定值
#include "zeroGradientFvPatchFields.H"  //边界条件，零梯度
#include "fixedFluxPressureFvPatchScalarField.H"  //边界条件，将压力梯度设置为所提供的值，从而使边界上的通量是由速度边界条件所规定的
#include "constrainHbyA.H"    //满足动量守恒，HbyA是方程求解过程中，移向相关的计算式子
#include "constrainPressure.H"
#include "adjustPhi.H"        //修正通量计算，调整通量的平衡以满足连续性。
#include "findRefCell.H"      //找到离给定单元最近的参考单元，貌似与probes有关
#include "IOMRFZoneList.H"    //MRF域列表读取
#include "constants.H"        //常数字典，pi、e等，算例配置constant文件相关内容，与物理场类型有关，如热、磁场、化学反应等
#include "gravityMeshObject.H"  //重力加速度，注册在Time实例里，对于计算域内所有网格均相同

#include "columnFvMesh.H"       //基于现有的网格和/或字段，生成网格的一维柱状表示

#include "OSspecific.H"       //与系统配置和类型有关
#include "argList.H"          //程序运行命令参数读取
#include "timeSelector.H"     //用于选择时间的scalarRange列表

#ifndef namespaceFoam
#define namespaceFoam
    using namespace Foam;
#endif

#endif
