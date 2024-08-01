#include <iostream>
#include <Eigen/Dense>
#include "data_bus.h"
#include "qpOASES.hpp"
// #include "Timer.h"

const uint16_t mpc_N = 2;              //MPC预测步长
const uint16_t mpc_U_N = 2;
const uint16_t nx = 12;
const uint16_t nu = 13;

//约束
const uint16_t ncfr = 4;                //约束力数量
// const uint16_t ncfr = 5;                //约束力数量
const uint16_t ncf = 4;                 //受力足数
const uint16_t nc = ncf * ncfr;


class MPC{
public:
    MPC(double dtIN);

    void set_weight(double u_weight, Eigen::MatrixXd L_diag, Eigen::MatrixXd K_diag);
    void cal();

    void dataBusRead(DataBus &Data);
    void dataBusWrite(DataBus &Data);

    void enable();
    void disable();
    bool get_ENA();

private:
    void copy_Eigen_to_real_t(qpOASES::real_t* target, Eigen::MatrixXd source, int nRows, int nCols);
    bool EN = false;

    Eigen::Matrix<double,nx,nx>   Ac[mpc_N], A[mpc_N];
    Eigen::Matrix<double,nx,nu>   Bc[mpc_N], B[mpc_N];
    Eigen::Matrix<double,nx,1>    Cc, C;

    Eigen::Matrix<double,nx*mpc_N,nx>         Aqp;
    Eigen::Matrix<double,nx*mpc_N,nx*mpc_N>   Aqp1;
    Eigen::Matrix<double,nx*mpc_N,nu*mpc_N>   Bqp1;
    Eigen::Matrix<double,nx*mpc_N,nu*mpc_U_N>      Bqp;
    Eigen::Matrix<double,nx*mpc_N,1>          Cqp1;
    Eigen::Matrix<double,nx*mpc_N,1>          Cqp;

    Eigen::Matrix<double,nu*mpc_U_N,1>           Ufe;
    Eigen::Matrix<double,nu,1>              Ufe_pre;
    Eigen::Matrix<double,nx*mpc_N,1>        Xd;                          //期望状态
    Eigen::Matrix<double,nx,1>              X_cur;                       //当前状态
    Eigen::Matrix<double,nx,1>              X_cal;
    Eigen::Matrix<double,nx,1>              X_cal_pre;
    Eigen::Matrix<double,nx,1>              dX_cal;

    Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic>    Q;
    Eigen::Matrix<double,nu*mpc_U_N, nu*mpc_U_N>            R, M;
    double alpha;
    Eigen::Matrix<double,nu*mpc_U_N, nu*mpc_U_N>        H;
    Eigen::Matrix<double,1, nu*mpc_U_N>                 c;
    Eigen::Matrix<double,nu*mpc_U_N,1>                  u_low, u_up;
    Eigen::Matrix<double,nc*mpc_U_N, nu*mpc_U_N>        As;
    Eigen::Matrix<double,nc*mpc_U_N,1>                  bs;
    double      max[6], min[6];

    double m, g, miu;
    Eigen::Matrix<double,3,1>   pCoM;
    Eigen::Matrix<double,12,1>  pf2com, pf2comd, pe;
    Eigen::Matrix<double,12,1>   pf2comi[mpc_N];    //机身到足端的向量表示的数组
    Eigen::Matrix<double,3,3>   Ic;                 //机器人质心的全身惯量矩阵
    Eigen::Matrix<double,3,3>   R_curz[mpc_N];
    Eigen::Matrix<double,3,3>   R_cur;
    Eigen::Matrix<double,3,3>   R_w2f, R_f2w;       //世界坐标系下左右脚旋转矩阵与转置

    int legStateCur;
    int legStateNext;
    int legState[10];
    double  dt;

    qpOASES::QProblem QP;
    qpOASES::real_t qp_H[nu*mpc_U_N * nu*mpc_U_N];   // H
    qpOASES::real_t qp_As[nc*mpc_U_N * nu*mpc_U_N];  // A（约束矩阵）
    qpOASES::real_t qp_c[nu*mpc_U_N];                // f
    qpOASES::real_t qp_lbA[nc*mpc_U_N];              // 约束下限
    qpOASES::real_t qp_ubA[nc*mpc_U_N];              // 约束上限
    qpOASES::real_t qp_lu[nu*mpc_U_N];               // 控制量下限
    qpOASES::real_t qp_uu[nu*mpc_U_N];               // 控制量上限
    qpOASES::int_t nWSR=100;                    // 最大重新计算次数
    qpOASES::real_t cpu_time=0.1;               // CPU限制时间
    qpOASES::real_t xOpt_iniGuess[nu*mpc_U_N];

	double			qp_cpuTime;
    int 			qp_Status, qp_nWSR;

};