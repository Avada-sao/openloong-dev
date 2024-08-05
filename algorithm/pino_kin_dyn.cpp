/*
This is part of OpenLoong Dynamics Control, an open project for the control of biped robot,
Copyright (C) 2024 Humanoid Robot (Shanghai) Co., Ltd, under Apache 2.0.
Feel free to use in any purpose, and cite OpenLoong-Dynamics-Control in any style, to contribute to the advancement of the community.
 <https://atomgit.com/openloong/openloong-dyn-control.git>
 <web@openloong.org.cn>
*/
#include "pino_kin_dyn.h"
#include <utility>

Pin_KinDyn::Pin_KinDyn(std::string urdf_pathIn) {
    pinocchio::JointModelFreeFlyer root_joint;//在buildModel()时调用pinocchio::JointModelFreeFlyer声明的root_joint，建立的模型是在urdf中第一个节点和世界坐标系之间加入了一个6自由度的浮动关节。
    pinocchio::urdf::buildModel(urdf_pathIn,root_joint,model_biped);
    pinocchio::urdf::buildModel(urdf_pathIn,model_biped_fixed);
    data_biped=pinocchio::Data(model_biped);
    data_biped_fixed=pinocchio::Data(model_biped_fixed);
    model_nv=model_biped.nv;

    std::cout<<"model name: "<<model_biped.name<<std::endl;
    for (pinocchio::JointIndex idx = 1; idx < model_biped.njoints; idx++)
        std::cout<<"Joint #"<<idx<<": "<<model_biped.names[idx]<<std::endl;

    J_fr=Eigen::MatrixXd::Zero(6,model_nv);
    J_fl=Eigen::MatrixXd::Zero(6,model_nv);
    J_rr=Eigen::MatrixXd::Zero(6,model_nv);
    J_rl=Eigen::MatrixXd::Zero(6,model_nv);
    J_base=Eigen::MatrixXd::Zero(6,model_nv);

    dJ_fr=Eigen::MatrixXd::Zero(6,model_nv);
    dJ_fl=Eigen::MatrixXd::Zero(6,model_nv);
    dJ_rr=Eigen::MatrixXd::Zero(6,model_nv);
    dJ_rl=Eigen::MatrixXd::Zero(6,model_nv);
    dJ_base=Eigen::MatrixXd::Zero(6,model_nv);

    q.setZero();
    dq.setZero();
    ddq.setZero();
    Rcur.setIdentity();
    dyn_M=Eigen::MatrixXd::Zero(model_nv,model_nv);
    dyn_M_inv=Eigen::MatrixXd::Zero(model_nv,model_nv);
    dyn_C=Eigen::MatrixXd::Zero(model_nv,model_nv);
    dyn_G=Eigen::MatrixXd::Zero(model_nv,1);

    // get joint index for Pinocchio Lib, need to redefined the joint name for new model
    fr_foot=model_biped.getFrameId("FR_foot");
    fl_foot=model_biped.getFrameId("FL_foot");
    rr_foot=model_biped.getFrameId("RR_foot");
    rl_foot=model_biped.getFrameId("RL_foot");

    fr_hip_joint=model_biped.getJointId("FR_hip_joint");
    fl_hip_joint=model_biped.getJointId("FL_hip_joint");
    rr_hip_joint=model_biped.getJointId("RR_hip_joint");
    rl_hip_joint=model_biped.getJointId("RL_hip_joint");
    base_joint=model_biped.getJointId("root_joint");

    // read joint pvt parameters
    //Json::Reader reader;
    //Json::Value root_read;
    //std::ifstream in("joint_ctrl_config.json",std::ios::binary);
    // motorMaxTorque=Eigen::VectorXd::Zero(motorName.size());
    // motorMaxPos=Eigen::VectorXd::Zero(motorName.size());
    // motorMinPos=Eigen::VectorXd::Zero(motorName.size());
    // reader.parse(in,root_read);
    // for (int i=0;i<motorName.size();i++){
    //     motorMaxTorque(i)=(root_read[motorName[i]]["maxTorque"].asDouble());
    //     motorMaxPos(i)=(root_read[motorName[i]]["maxPos"].asDouble());
    //     motorMinPos(i)=(root_read[motorName[i]]["minPos"].asDouble());
    // }

    motorReachLimit.assign(motorName.size(),false);
    tauJointOld=Eigen::VectorXd::Zero(motorName.size());
}

void Pin_KinDyn::dataBusRead(const DataBus &robotState) {
    //  For Pinocchio: The base translation part is expressed in the parent frame (here the world coordinate system)
    //  while its velocity is expressed in the body coordinate system.
    //  https://github.com/stack-of-tasks/pinocchio/issues/1137
    //  q = [global_base_position, global_base_quaternion, joint_positions]
    //  v = [local_base_velocity_linear, local_base_velocity_angular, joint_velocities]
    q=robotState.q;
    dq=robotState.dq;
    dq.block(0,0,3,1)=robotState.base_rot.transpose()*dq.block(0,0,3,1);
    dq.block(3,0,3,1)=robotState.base_rot.transpose()*dq.block(3,0,3,1);
    ddq=robotState.ddq;
}

void Pin_KinDyn::dataBusWrite(DataBus &robotState) {
    // NOTE: for the following Jacobians, x=J*dq, both x and dq are defined in the world frame
    robotState.J_fr=J_fr;
    robotState.J_fl=J_fl;
    robotState.J_rr=J_rr;
    robotState.J_rl=J_rl;
    robotState.dJ_fr=dJ_fr;
    robotState.dJ_fl=dJ_fl;
    robotState.dJ_rr=dJ_rr;
    robotState.dJ_rl=dJ_rl;

    robotState.J_base=J_base;
    robotState.dJ_base=dJ_base;

    robotState.fe_fr_pos_W=fe_fr_pos;
    robotState.fe_fl_pos_W=fe_fl_pos;
    robotState.fe_rr_pos_W=fe_rr_pos;
    robotState.fe_rl_pos_W=fe_rl_pos;

    robotState.fe_l_pos_L=fe_l_pos_body;
    robotState.fe_r_pos_L=fe_r_pos_body;

    robotState.hip_r_pos_L=hip_r_pos_body;
    robotState.hip_l_pos_L=hip_l_pos_body;
    robotState.hip_fr_pos_W=hip_fr_pos;
    robotState.hip_fl_pos_W=hip_fl_pos;
    robotState.hip_rr_pos_W=hip_rr_pos;
    robotState.hip_rl_pos_W=hip_rl_pos;

    robotState.dyn_M=dyn_M;
    robotState.dyn_M_inv=dyn_M_inv;
    robotState.dyn_C=dyn_C;
    robotState.dyn_G=dyn_G;
    robotState.dyn_Ag=dyn_Ag;
    robotState.dyn_dAg=dyn_dAg;
    robotState.dyn_Non=dyn_Non;

    robotState.inertia = inertia;  // w.r.t body frame
}

// update jacobians and joint positions
void Pin_KinDyn::computeJ_dJ() {
    
    pinocchio::forwardKinematics(model_biped,data_biped,q);
    pinocchio::jacobianCenterOfMass(model_biped, data_biped, q, true);//在世界坐标系下质心位置的雅各布矩阵雅可比矩阵，(matrix 3 x model.nv)
    pinocchio::computeJointJacobians(model_biped,data_biped,q);
    pinocchio::computeJointJacobiansTimeVariation(model_biped,data_biped,q,dq);//雅各布随时间的变化，dj (matrix 6 x model.nv)
    pinocchio::updateGlobalPlacements(model_biped,data_biped);//更新关节位置

    pinocchio::getFrameJacobian(model_biped,data_biped,fr_foot,pinocchio::LOCAL_WORLD_ALIGNED,J_fr);//LOCAL_WORLD_ALIGNED:坐标原点与机身坐标系重合，但坐标系轴线与世界坐标系平行
    pinocchio::getFrameJacobian(model_biped,data_biped,fl_foot,pinocchio::LOCAL_WORLD_ALIGNED,J_fl);//(matrix 6 x model.nv)
    pinocchio::getFrameJacobian(model_biped,data_biped,rr_foot,pinocchio::LOCAL_WORLD_ALIGNED,J_rr);
    pinocchio::getFrameJacobian(model_biped,data_biped,rl_foot,pinocchio::LOCAL_WORLD_ALIGNED,J_rl);
    pinocchio::getJointJacobian(model_biped,data_biped,base_joint,pinocchio::LOCAL_WORLD_ALIGNED,J_base);

    pinocchio::getFrameJacobianTimeVariation(model_biped,data_biped,fr_foot,pinocchio::LOCAL_WORLD_ALIGNED,dJ_fr);
    pinocchio::getFrameJacobianTimeVariation(model_biped,data_biped,fl_foot,pinocchio::LOCAL_WORLD_ALIGNED,dJ_fl);
    pinocchio::getFrameJacobianTimeVariation(model_biped,data_biped,rr_foot,pinocchio::LOCAL_WORLD_ALIGNED,dJ_rr);
    pinocchio::getFrameJacobianTimeVariation(model_biped,data_biped,rl_foot,pinocchio::LOCAL_WORLD_ALIGNED,dJ_rl);
    pinocchio::getJointJacobianTimeVariation(model_biped,data_biped,base_joint,pinocchio::LOCAL_WORLD_ALIGNED,dJ_base);

    // fe_fr_pos=data_biped.oMi[fr_foot].translation();//oMi:物体相对于参考坐标系的运动变换
    // fe_fl_pos=data_biped.oMi[fl_foot].translation(); 
    // fe_rr_pos=data_biped.oMi[rr_foot].translation();
    // fe_rl_pos=data_biped.oMi[rl_foot].translation();  
    hip_fl_pos=data_biped.oMi[fl_hip_joint].translation();
    hip_fr_pos=data_biped.oMi[fr_hip_joint].translation();
    hip_rl_pos=data_biped.oMi[rl_hip_joint].translation();
    hip_rr_pos=data_biped.oMi[rr_hip_joint].translation();
    base_pos=data_biped.oMi[base_joint].translation();
    base_rot=data_biped.oMi[base_joint].rotation();

    legPos=Eigen::MatrixXd::Zero(3,4);
    for(int i=0; i<4; i++){
        Eigen::Vector3d legJoint;
        legJoint = q.block(7+i*3,0,3,1);
        legPos.block(0,i,3,1) = computeFp(legJoint, i);
    }

    fe_fl_pos = hip_fl_pos + legPos.block(0,0,3,1);
    fe_fr_pos = hip_fr_pos + legPos.block(0,1,3,1);
    fe_rl_pos = hip_rl_pos + legPos.block(0,2,3,1);
    fe_rr_pos = hip_rr_pos + legPos.block(0,3,3,1);

    // std::cout<<"base_pos"<<base_pos.transpose()<<std::endl;
    // std::cout<<"fe_fr_pos"<<fe_fr_pos.transpose()<<std::endl;
    // std::cout<<"fe_fl_pos"<<fe_fl_pos.transpose()<<std::endl;
    // std::cout<<"fe_rr_pos"<<fe_rr_pos.transpose()<<std::endl;
    // std::cout<<"fe_rl_pos"<<fe_rl_pos.transpose()<<std::endl;

    Jcom=data_biped.Jcom;// v_{\text{CoM}} = J_{\text{CoM}} \dot{q}

    Eigen::MatrixXd Mpj; // transform into world frame, and accept dq that in world frame
    Mpj=Eigen::MatrixXd::Identity(model_nv,model_nv);
    Mpj.block(0,0,3,3)=base_rot.transpose();
    Mpj.block(3,3,3,3)=base_rot.transpose();
    J_fr=J_fr*Mpj;
    J_fl=J_fl*Mpj;
    J_rr=J_rr*Mpj;
    J_rl=J_rl*Mpj;   
    J_base=J_base*Mpj;

    dJ_fr=dJ_fr*Mpj;
    dJ_fl=dJ_fl*Mpj;
    dJ_rr=dJ_rr*Mpj;
    dJ_rl=dJ_rl*Mpj;
    dJ_base=dJ_base*Mpj;
    Jcom=Jcom*Mpj;

    // Eigen::VectorXd q_fixed;
    // q_fixed=q.block(7,0,model_biped_fixed.nv,1);
    // pinocchio::forwardKinematics(model_biped_fixed,data_biped_fixed,q_fixed);
    // pinocchio::updateGlobalPlacements(model_biped_fixed,data_biped_fixed);
    // fe_l_pos_body=data_biped_fixed.oMi[l_ankle_joint_fixed].translation();
    // fe_r_pos_body=data_biped_fixed.oMi[r_ankle_joint_fixed].translation();
    // hip_l_pos_body=data_biped_fixed.oMi[l_hip_joint_fixed].translation();
    // hip_r_pos_body=data_biped_fixed.oMi[r_hip_joint_fixed].translation();
}

Eigen::Vector3d Pin_KinDyn::computeFp(Eigen::Vector3d q, int leg) {

    double sideSigns[4] = {1, -1, 1, -1};
    double sideSign = sideSigns[leg];

    double l1 = 0.08505;
    double l2 = 0.2;
    double l3 = 0.2;

    double s1 = std::sin(q(0));
    double s2 = std::sin(q(1));
    double s3 = std::sin(q(2));

    double c1 = std::cos(q(0));
    double c2 = std::cos(q(1));
    double c3 = std::cos(q(2));

    double c23 = c2 * c3 - s2 * s3;
    double s23 = s2 * c3 + c2 * s3;

    Eigen::Vector3d p;
    p[0] = l3 * s23 + l2 * s2;
    p[1] = l1 * sideSign * c1 + l3 * (s1 * c23) + l2 * c2 * s1;
    p[2] = l1 * sideSign * s1 - l3 * (c1 * c23) - l2 * c1 * c2;
    return p;
}


Eigen::Quaterniond Pin_KinDyn::intQuat(const Eigen::Quaterniond &quat, const Eigen::Matrix<double, 3, 1> &w) {
    Eigen::Matrix3d Rcur=quat.normalized().toRotationMatrix();
    Eigen::Matrix3d Rinc=Eigen::Matrix3d::Identity();
    double theta=w.norm();
    if (theta>1e-8) {
        Eigen::Vector3d w_norm;
        w_norm = w / theta;
        Eigen::Matrix3d a;
        a << 0, -w_norm(2), w_norm(1),
                w_norm(0), 0, -w_norm(0),
                -w_norm(1), w_norm(0), 0;
        Rinc=Eigen::Matrix3d::Identity()+a*sin(theta)+a*a*(1-cos(theta));
    }
    Eigen::Matrix3d Rend=Rcur*Rinc;
    Eigen::Quaterniond quatRes;
    quatRes=Rend;
    return quatRes;
}

// intergrate the q with dq, for floating base dynamics
Eigen::VectorXd Pin_KinDyn::integrateDIY(const Eigen::VectorXd &qI, const Eigen::VectorXd &dqI) {
    Eigen::VectorXd qRes=Eigen::VectorXd::Zero(model_nv+1);
    Eigen::Vector3d wDes;
    wDes<<dqI(3),dqI(4),dqI(5);
    Eigen::Quaterniond quatNew,quatNow;
    quatNow.x()=qI(3);
    quatNow.y()=qI(4);
    quatNow.z()=qI(5);
    quatNow.w()=qI(6);
    quatNew= intQuat(quatNow,wDes);
    qRes=qI;
    qRes(0)+=dqI(0);qRes(1)+=dqI(1);qRes(2)+=dqI(2);
    qRes(3)=quatNew.x();qRes(4)=quatNew.y();qRes(5)=quatNew.z();qRes(6)=quatNew.w();
    for (int i=0;i<model_nv-6;i++){
        // if (i==0 || i==3){
        //     qRes(7) += dqI(9);
        //     qRes(10) += dqI(6);
        // }
        qRes(7+i)+=dqI(6+i);
    }

    

    return qRes;
}

// update dynamic parameters, M*ddq+C*dq+G=tau
void Pin_KinDyn::computeDyn() {
    // cal M
    pinocchio::crba(model_biped, data_biped, q);//关节惯性矩阵M上三角 部分（model.nv*model.nv）
    // Pinocchio only gives half of the M, needs to restore it here
    data_biped.M.triangularView<Eigen::Lower>() = data_biped.M.transpose().triangularView<Eigen::Lower>();
    dyn_M = data_biped.M;

    // cal Minv
    pinocchio::computeMinverse(model_biped, data_biped, q);
    data_biped.Minv.triangularView<Eigen::Lower>() = data_biped.Minv.transpose().triangularView<Eigen::Lower>();
    dyn_M_inv=data_biped.Minv;

    // cal C
    pinocchio::computeCoriolisMatrix(model_biped, data_biped, q, dq);
    dyn_C = data_biped.C;

    // cal G
    pinocchio::computeGeneralizedGravity(model_biped, data_biped, q);
    dyn_G = data_biped.g;

    // cal Ag, Centroidal Momentum Matrix. First three rows: linear, other three rows: angular
    pinocchio::dccrba(model_biped,data_biped,q,dq);//dM
    pinocchio::computeCentroidalMomentum(model_biped,data_biped,q,dq);
    dyn_Ag=data_biped.Ag;
    dyn_dAg=data_biped.dAg;

    // cal nonlinear item
    dyn_Non=dyn_C*dq+dyn_G;

    // cal I
    pinocchio::ccrba(model_biped, data_biped, q, dq);
    inertia = data_biped.Ig.inertia().matrix();//inertia:mpc中的简化模型的转动惯量矩阵（3*3）

    // cal CoM
    CoM_pos = data_biped.com[0];
    // std::cout<<"CoM_W"<<std::endl;
    // std::cout<<CoM_pos.transpose()<<std::endl;

    Eigen::MatrixXd Mpj, Mpj_inv; // transform into world frame
    Mpj=Eigen::MatrixXd::Identity(model_nv,model_nv);
    Mpj_inv=Eigen::MatrixXd::Identity(model_nv,model_nv);
    Mpj.block(0,0,3,3)=base_rot.transpose();
    Mpj.block(3,3,3,3)=base_rot.transpose();
    Mpj_inv.block(0,0,3,3)=base_rot;
    Mpj_inv.block(3,3,3,3)=base_rot;
    dyn_M=Mpj_inv*dyn_M*Mpj;
    dyn_M_inv=Mpj_inv*dyn_M_inv*Mpj;
    dyn_C=Mpj_inv*dyn_C*Mpj;
    dyn_G=Mpj_inv*dyn_G;
    dyn_Non=Mpj_inv*dyn_Non;
}

// Inverse kinematics for leg posture. Note: the Rdes and Pdes are both w.r.t the baselink coordinate in body frame!
Pin_KinDyn::IkRes
Pin_KinDyn::computeInK_Leg(const Eigen::Matrix3d &Rdes_L, const Eigen::Vector3d &Pdes_L, const Eigen::Matrix3d &Rdes_R,
                           const Eigen::Vector3d &Pdes_R) {
    const pinocchio::SE3 oMdesL(Rdes_L, Pdes_L);
    const pinocchio::SE3 oMdesR(Rdes_R, Pdes_R);
    // arm-l: 0-6, arm-r: 7-13, head: 14,15 waist: 16-18, leg-l: 19-24, leg-r: 25-30
    Eigen::VectorXd qIk=Eigen::VectorXd::Zero(model_biped_fixed.nv); // initial guess
    qIk[22]=-0.1;
    qIk[28]=-0.1;

    const double eps  = 1e-4;
    const int IT_MAX  = 100;
    const double DT   = 7e-1;
    const double damp = 5e-3;
    Eigen::MatrixXd JL(6,model_biped_fixed.nv);
    Eigen::MatrixXd JR(6,model_biped_fixed.nv);
    Eigen::MatrixXd JCompact(12,model_biped_fixed.nv);
    JL.setZero();
    JR.setZero();
    JCompact.setZero();

    bool success = false;
    Eigen::Matrix<double, 6, 1> errL, errR;
    Eigen::Matrix<double, 12,1> errCompact;
    Eigen::VectorXd v(model_biped_fixed.nv);

    pinocchio::JointIndex J_Idx_l, J_Idx_r;
    // J_Idx_l = l_ankle_joint_fixed;
    // J_Idx_r = r_ankle_joint_fixed;
    int itr_count{0};
    for (itr_count=0;; itr_count++)
    {
        pinocchio::forwardKinematics(model_biped_fixed,data_biped_fixed,qIk);
        const pinocchio::SE3 iMdL = data_biped_fixed.oMi[J_Idx_l].actInv(oMdesL);
        const pinocchio::SE3 iMdR = data_biped_fixed.oMi[J_Idx_r].actInv(oMdesR);
        errL = pinocchio::log6(iMdL).toVector();  // in joint frame
        errR = pinocchio::log6(iMdR).toVector();  // in joint frame
        errCompact.block<6,1>(0,0)=errL;
        errCompact.block<6,1>(6,0)=errR;
        if(errCompact.norm() < eps)
        {
            success = true;
            break;
        }
        if (itr_count >= IT_MAX)
        {
            success = false;
            break;
        }

        pinocchio::computeJointJacobian(model_biped_fixed,data_biped_fixed,qIk,J_Idx_l,JL);  // JL in joint frame
        pinocchio::computeJointJacobian(model_biped_fixed,data_biped_fixed,qIk,J_Idx_r,JR);  // JR in joint frame
        Eigen::MatrixXd W;
        W=Eigen::MatrixXd::Identity(model_biped_fixed.nv, model_biped_fixed.nv);   // weighted matrix
        // arm-l: 0-6, arm-r: 7-13, head: 14,15 waist: 16-18, leg-l: 19-24, leg-r: 25-30
        // W(16,16)=0.001;  // use a smaller value to make the solver try not to use waist joint
        // W(17,17)=0.001;
        // W(18,18)=0.001;
        JL.block(0,16,6,3).setZero();
        JR.block(0,16,6,3).setZero();
        pinocchio::Data::Matrix6 JlogL;
        pinocchio::Data::Matrix6 JlogR;
        pinocchio::Jlog6(iMdL.inverse(), JlogL);
        pinocchio::Jlog6(iMdR.inverse(), JlogR);
        JL = -JlogL * JL;
        JR = -JlogR * JR;
        JCompact.block(0,0,6,model_biped_fixed.nv)=JL;
        JCompact.block(6,0,6,model_biped_fixed.nv)=JR;
        // pinocchio::Data::Matrix6 JJt;
        Eigen::Matrix<double,12,12> JJt;
        JJt.noalias() = JCompact * W * JCompact.transpose();
        JJt.diagonal().array() += damp;
        v.noalias() = - W * JCompact.transpose() * JJt.ldlt().solve(errCompact);
        qIk = pinocchio::integrate(model_biped_fixed,qIk,v*DT);
    }

    IkRes res;
    res.err=errCompact;
    res.itr=itr_count;

    if(success){
        res.status=0;
    }
    else{
        res.status=-1;
    }
    res.jointPosRes=qIk;
    return res;
}

// Inverse Kinematics for hand posture. Note: the Rdes and Pdes are both w.r.t the baselink coordinate in body frame!
Pin_KinDyn::IkRes
Pin_KinDyn::computeInK_Hand(const Eigen::Matrix3d &Rdes_L, const Eigen::Vector3d &Pdes_L, const Eigen::Matrix3d &Rdes_R,
                            const Eigen::Vector3d &Pdes_R) {
    const pinocchio::SE3 oMdesL(Rdes_L, Pdes_L);
    const pinocchio::SE3 oMdesR(Rdes_R, Pdes_R);
    Eigen::VectorXd qIk=Eigen::VectorXd::Zero(model_biped_fixed.nv); // initial guess
    // arm-l: 0-6, arm-r: 7-13, head: 14,15 waist: 16-18, leg-l: 19-24, leg-r: 25-30
    qIk.block<7,1>(0,0)<< 0.433153883479341,    1.11739345867607,    1.88491913406236,
            0.802378252758275,    1.22726400279662,    0.0249797771339966,  -0.0875282610654057;

    qIk.block<7,1>(7,0)<<-0.433152540054138,   -1.11739347975224,  -1.88492038240761,
            0.802375980602373,   -1.22726323451626,   0.0249795712262396, 0.0875271396314979;

    const double eps  = 1e-4;
    const int IT_MAX  = 100;
    const double DT   = 6e-1;
    const double damp = 1e-2;
    Eigen::MatrixXd JL(6,model_biped_fixed.nv);
    Eigen::MatrixXd JR(6,model_biped_fixed.nv);
    Eigen::MatrixXd JCompact(12,model_biped_fixed.nv);
    JL.setZero();
    JR.setZero();
    JCompact.setZero();

    bool success = false;
    Eigen::Matrix<double, 6, 1> errL, errR;
    Eigen::Matrix<double, 12,1> errCompact;
    Eigen::VectorXd v(model_biped_fixed.nv);

    pinocchio::JointIndex J_Idx_l, J_Idx_r;
    // J_Idx_l = l_hand_joint_fixed;
    // J_Idx_r = r_hand_joint_fixed;
    int itr_count{0};
    for (itr_count=0;; itr_count++)
    {
        pinocchio::forwardKinematics(model_biped_fixed,data_biped_fixed,qIk);
        const pinocchio::SE3 iMdL = data_biped_fixed.oMi[J_Idx_l].actInv(oMdesL);
        const pinocchio::SE3 iMdR = data_biped_fixed.oMi[J_Idx_r].actInv(oMdesR);
        errL = pinocchio::log6(iMdL).toVector();  // in joint frame
        errR = pinocchio::log6(iMdR).toVector();  // in joint frame
        errCompact.block<6,1>(0,0)=errL;
        errCompact.block<6,1>(6,0)=errR;
        if(errCompact.norm() < eps)
        {
            success = true;
            break;
        }
        if (itr_count >= IT_MAX)
        {
            success = false;
            break;
        }

        pinocchio::computeJointJacobian(model_biped_fixed,data_biped_fixed,qIk,J_Idx_l,JL);  // JL in joint frame
        pinocchio::computeJointJacobian(model_biped_fixed,data_biped_fixed,qIk,J_Idx_r,JR);  // JR in joint frame
        pinocchio::Data::Matrix6 JlogL;
        pinocchio::Data::Matrix6 JlogR;
        pinocchio::Jlog6(iMdL.inverse(), JlogL);
        pinocchio::Jlog6(iMdR.inverse(), JlogR);
        JL = -JlogL * JL;
        JR = -JlogR * JR;
        JCompact.block(0,0,6,model_biped_fixed.nv)=JL;
        JCompact.block(6,0,6,model_biped_fixed.nv)=JR;
        // pinocchio::Data::Matrix6 JJt;
        Eigen::Matrix<double,12,12> JJt;
        JJt.noalias() = JCompact * JCompact.transpose();
        JJt.diagonal().array() += damp;
        v.noalias() = - JCompact.transpose() * JJt.ldlt().solve(errCompact);
        qIk = pinocchio::integrate(model_biped_fixed,qIk,v*DT);
    }

    IkRes res;
    res.err=errCompact;
    res.itr=itr_count;

    if(success){
        res.status=0;
    }
    else{
        res.status=-1;
    }

    res.jointPosRes=qIk;
    return res;
}

// must call computeDyn() first!
void Pin_KinDyn::workspaceConstraint(Eigen::VectorXd &qFT, Eigen::VectorXd &tauJointFT) {
    for (int i=0; i<motorName.size();i++)
        if (qFT(i+7)>motorMaxPos(i)){
            qFT(i+7)=motorMaxPos(i);
            motorReachLimit[i]=true;
            tauJointFT(i)=tauJointOld(i);
        }
        else if (qFT(i+7)<motorMinPos(i)){
            qFT(i+7)=motorMinPos(i);
            motorReachLimit[i]=true;
            tauJointFT(i)=tauJointOld(i);
        }
        else
            motorReachLimit[i]=false;
    tauJointOld=tauJointFT;
}

























