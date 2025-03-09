
#pragma once

#include <array>

#include <stdexcept>

template<typename T>
struct Spacetime { 

   Spacetime() :
      alpha_(U[0]), 
      beta_X(U[1]), beta_Y(U[2]), beta_Z(U[3]), 
      gamma_xx(U[4]), gamma_xy(U[5]), gamma_xz(U[6]), 
      gamma_yy(U[7]), gamma_yz(U[8]), 
      gamma_zz(U[9]), 
      K_xx(U[10]), K_xy(U[11]), K_xz(U[12]), 
      K_yy(U[13]), K_yz(U[14]), 
      K_zz(U[15]), 
      A_x(U[16]), A_y(U[17]), A_z(U[18]), 
      B_xX(U[19]), B_xY(U[20]), B_xZ(U[21]), 
      B_yX(U[22]), B_yY(U[23]), B_yZ(U[24]), 
      B_zX(U[25]), B_zY(U[26]), B_zZ(U[27]), 
      D_xxx(U[28]), D_xxy(U[29]), D_xxz(U[30]), 
      D_xyy(U[31]), D_xyz(U[32]), 
      D_xzz(U[33]), 
      D_yxx(U[34]), D_yxy(U[35]), D_yxz(U[36]), 
      D_yyy(U[37]), D_yyz(U[38]), 
      D_yzz(U[39]), 
      D_zxx(U[40]), D_zxy(U[41]), D_zxz(U[42]), 
      D_zyy(U[43]), D_zyz(U[44]), 
      D_zzz(U[45]), 
      Theta_(U[46]), 
      Z_x(U[47]), Z_y(U[48]), Z_z(U[49])
   {}
   
   Spacetime flux(int dim) const {
      Spacetime F;
      
      T gamma_XX = gamma_yy * gamma_zz - gamma_yz * gamma_yz;
      T gamma_YY = gamma_xx * gamma_zz - gamma_xz * gamma_xz;
      T gamma_ZZ = gamma_xx * gamma_yy - gamma_xy * gamma_xy;
      T gamma_XY = gamma_zz * gamma_xy - gamma_xz * gamma_yz;
      T gamma_XZ = gamma_xy * gamma_yz - gamma_yy * gamma_xz;
      T gamma_YZ = gamma_yz * gamma_xx - gamma_xz * gamma_xy;
      
      T const  gamma_ = gamma_xx * gamma_XX + gamma_xy * gamma_XY + gamma_xz * gamma_XZ;
      
      T const  igamma_ = T(1) / gamma_;
      
      gamma_XX *= igamma_;
      gamma_YY *= igamma_;
      gamma_ZZ *= igamma_;
      gamma_XY *= igamma_;
      gamma_XZ *= igamma_;
      gamma_YZ *= igamma_;
      
      T const  D_x = D_xxx * gamma_XX + D_xyy * gamma_YY + D_xzz * gamma_ZZ + T(2) * (D_xxy * gamma_XY + D_xxz * gamma_XZ + D_xyz * gamma_YZ);
      T const  D_y = D_yxx * gamma_XX + D_yyy * gamma_YY + D_yzz * gamma_ZZ + T(2) * (D_yxy * gamma_XY + D_yxz * gamma_XZ + D_yyz * gamma_YZ);
      T const  D_z = D_zxx * gamma_XX + D_zyy * gamma_YY + D_zzz * gamma_ZZ + T(2) * (D_zxy * gamma_XY + D_zxz * gamma_XZ + D_zyz * gamma_YZ);
      
      T const  E_x = D_xxx * gamma_XX + D_xxy * gamma_XY + D_xxz * gamma_XZ + D_yxx * gamma_XY + D_yxy * gamma_YY + D_yxz * gamma_YZ + D_zxx * gamma_XZ + D_zxy * gamma_YZ + D_zxz * gamma_ZZ;
      T const  E_y = D_xxy * gamma_XX + D_xyy * gamma_XY + D_xyz * gamma_XZ + D_yxy * gamma_XY + D_yyy * gamma_YY + D_yyz * gamma_YZ + D_zxy * gamma_XZ + D_zyy * gamma_YZ + D_zyz * gamma_ZZ;
      T const  E_z = D_xxz * gamma_XX + D_xyz * gamma_XY + D_xzz * gamma_XZ + D_yxz * gamma_XY + D_yyz * gamma_YY + D_yzz * gamma_YZ + D_zxz * gamma_XZ + D_zyz * gamma_YZ + D_zzz * gamma_ZZ;
      
      T const  V_x = D_x - E_x - Z_x;
      T const  V_y = D_y - E_y - Z_y;
      T const  V_z = D_z - E_z - Z_z;
      
      T const  K_ = gamma_XX * K_xx + gamma_YY * K_yy + gamma_ZZ * K_zz + T(2) * (gamma_XY * K_xy + gamma_XZ * K_xz + gamma_YZ * K_yz);
      
      T const  ahpla_ = T(1) / alpha_;
      
      T const  Q_ = K_ - T(2) * Theta_;
      
      T const  Q_x = alpha_ * (A_x - D_x + T(2) * V_x);
      T const  Q_y = alpha_ * (A_y - D_y + T(2) * V_y);
      T const  Q_z = alpha_ * (A_z - D_z + T(2) * V_z);
      
      T const  Q_xx = K_xx - ahpla_ * B_xX;
      T const  Q_yy = K_yy - ahpla_ * B_yY;
      T const  Q_zz = K_zz - ahpla_ * B_zZ;
      T const  Q_xy = K_xy - T(0.5) * ahpla_ * (B_xY + B_yX);
      T const  Q_xz = K_xz - T(0.5) * ahpla_ * (B_xZ + B_zX);
      T const  Q_yz = K_yz - T(0.5) * ahpla_ * (B_yZ + B_zY);
      
      T const  K_Xx = gamma_XX * K_xx + gamma_XY * K_xy + gamma_XZ * K_xz;
      T const  K_Xy = gamma_XX * K_xy + gamma_XY * K_yy + gamma_XZ * K_yz;
      T const  K_Xz = gamma_XX * K_xz + gamma_XY * K_yz + gamma_XZ * K_zz;
      T const  K_Yx = gamma_XY * K_xx + gamma_YY * K_xy + gamma_YZ * K_xz;
      T const  K_Yy = gamma_XY * K_xy + gamma_YY * K_yy + gamma_YZ * K_yz;
      T const  K_Yz = gamma_XY * K_xz + gamma_YY * K_yz + gamma_YZ * K_zz;
      T const  K_Zx = gamma_XZ * K_xx + gamma_YZ * K_xy + gamma_ZZ * K_xz;
      T const  K_Zy = gamma_XZ * K_xy + gamma_YZ * K_yy + gamma_ZZ * K_yz;
      T const  K_Zz = gamma_XZ * K_xz + gamma_YZ * K_yz + gamma_ZZ * K_zz;
      
      switch( dim ) {
      
      case 0: {
         
         T const  D_xXX = gamma_XX * (D_xxx * gamma_XX + D_xxy * gamma_XY + D_xxz * gamma_XZ) + gamma_XY * (D_xxy * gamma_XX + D_xyy * gamma_XY + D_xyz * gamma_XZ) + gamma_XZ * (D_xxz * gamma_XX + D_xyz * gamma_XY + D_xzz * gamma_XZ);
         T const  D_xXY = gamma_XX * (D_xxx * gamma_XY + D_xxy * gamma_YY + D_xxz * gamma_YZ) + gamma_XY * (D_xxy * gamma_XY + D_xyy * gamma_YY + D_xyz * gamma_YZ) + gamma_XZ * (D_xxz * gamma_XY + D_xyz * gamma_YY + D_xzz * gamma_YZ);
         T const  D_xXZ = gamma_XX * (D_xxx * gamma_XZ + D_xxy * gamma_YZ + D_xxz * gamma_ZZ) + gamma_XY * (D_xxy * gamma_XZ + D_xyy * gamma_YZ + D_xyz * gamma_ZZ) + gamma_XZ * (D_xxz * gamma_XZ + D_xyz * gamma_YZ + D_xzz * gamma_ZZ);
         T const  D_xYY = gamma_XY * (D_xxx * gamma_XY + D_xxy * gamma_YY + D_xxz * gamma_YZ) + gamma_YY * (D_xxy * gamma_XY + D_xyy * gamma_YY + D_xyz * gamma_YZ) + gamma_YZ * (D_xxz * gamma_XY + D_xyz * gamma_YY + D_xzz * gamma_YZ);
         T const  D_xYZ = gamma_XY * (D_xxx * gamma_XZ + D_xxy * gamma_YZ + D_xxz * gamma_ZZ) + gamma_YY * (D_xxy * gamma_XZ + D_xyy * gamma_YZ + D_xyz * gamma_ZZ) + gamma_YZ * (D_xxz * gamma_XZ + D_xyz * gamma_YZ + D_xzz * gamma_ZZ);
         T const  D_xZZ = gamma_XZ * (D_xxx * gamma_XZ + D_xxy * gamma_YZ + D_xxz * gamma_ZZ) + gamma_YZ * (D_xxy * gamma_XZ + D_xyy * gamma_YZ + D_xyz * gamma_ZZ) + gamma_ZZ * (D_xxz * gamma_XZ + D_xyz * gamma_YZ + D_xzz * gamma_ZZ);
         
         T const  lambda_xxx = D_xxx + A_x + D_x - T(2) * (E_x + Z_x);
         T const  lambda_xxy = D_xxy + T(0.5) * (A_y + D_y - T(2) * (E_y + Z_y));
         T const  lambda_xxz = D_xxz + T(0.5) * (A_z + D_z - T(2) * (E_z + Z_z));
         T const &lambda_xyy = D_xyy;
         T const &lambda_xyz = D_xyz;
         T const &lambda_xzz = D_xzz;
         T const &lambda_yxx = D_yxx;
         T const  lambda_yxy = D_yxy + T(0.5) * (A_x + D_x - T(2) * (E_x + Z_x));
         T const &lambda_yxz = D_yxz;
         T const  lambda_yyy = D_yyy + A_y + D_y - T(2) * (E_y + Z_y);
         T const  lambda_yyz = D_yyz + T(0.5) * (A_z + D_z - T(2) * (E_z + Z_z));
         T const &lambda_yzz = D_yzz;
         T const &lambda_zxx = D_zxx;
         T const &lambda_zxy = D_zxy;
         T const  lambda_zxz = D_zxz + T(0.5) * (A_x + D_x - T(2) * (E_x + Z_x));
         T const &lambda_zyy = D_zyy;
         T const  lambda_zyz = D_zyz + T(0.5) * (A_y + D_y - T(2) * (E_y + Z_y));
         T const  lambda_zzz = D_zzz + A_z + D_z - T(2) * (E_z + Z_z);
         T const  lambda_Xxx = gamma_XX * lambda_xxx + gamma_XY * lambda_yxx + gamma_XZ * lambda_zxx;
         T const  lambda_Xxy = gamma_XX * lambda_xxy + gamma_XY * lambda_yxy + gamma_XZ * lambda_zxy;
         T const  lambda_Xxz = gamma_XX * lambda_xxz + gamma_XY * lambda_yxz + gamma_XZ * lambda_zxz;
         T const  lambda_Xyy = gamma_XX * lambda_xyy + gamma_XY * lambda_yyy + gamma_XZ * lambda_zyy;
         T const  lambda_Xyz = gamma_XX * lambda_xyz + gamma_XY * lambda_yyz + gamma_XZ * lambda_zyz;
         T const  lambda_Xzz = gamma_XX * lambda_xzz + gamma_XY * lambda_yzz + gamma_XZ * lambda_zzz;
         
         F.alpha_ = F.beta_X = F.beta_Y = F.beta_Z = F.gamma_xx = F.gamma_xy = F.gamma_xz = F.gamma_yy = F.gamma_yz = F.gamma_zz = T(0);
         
         F.K_xx = -beta_X * K_xx + alpha_ * lambda_Xxx;
         F.K_xy = -beta_X * K_xy + alpha_ * lambda_Xxy;
         F.K_xz = -beta_X * K_xz + alpha_ * lambda_Xxz;
         F.K_yy = -beta_X * K_yy + alpha_ * lambda_Xyy;
         F.K_yz = -beta_X * K_yz + alpha_ * lambda_Xyz;
         F.K_zz = -beta_X * K_zz + alpha_ * lambda_Xzz;
         
         F.A_x = -beta_X * A_x + alpha_ * Q_;
         F.A_y = -beta_X * A_y;
         F.A_z = -beta_X * A_z;
         
         F.B_xX = -beta_X * B_xX + alpha_ * (gamma_XX * Q_x + gamma_XY * Q_y + gamma_XZ * Q_z);
         F.B_xY = -beta_X * B_xY + alpha_ * (gamma_XY * Q_x + gamma_YY * Q_y + gamma_YZ * Q_z);
         F.B_xZ = -beta_X * B_xZ + alpha_ * (gamma_XZ * Q_x + gamma_YZ * Q_y + gamma_ZZ * Q_z);
         F.B_yX = -beta_X * B_yX;
         F.B_yY = -beta_X * B_yY;
         F.B_yZ = -beta_X * B_yZ;
         F.B_zX = -beta_X * B_zX;
         F.B_zY = -beta_X * B_zY;
         F.B_zZ = -beta_X * B_zZ;
         
         F.D_xxx = -beta_X * D_xxx + alpha_ * Q_xx;
         F.D_xxy = -beta_X * D_xxy + alpha_ * Q_xy;
         F.D_xxz = -beta_X * D_xxz + alpha_ * Q_xz;
         F.D_xyy = -beta_X * D_xyy + alpha_ * Q_yy;
         F.D_xyz = -beta_X * D_xyz + alpha_ * Q_yz;
         F.D_xzz = -beta_X * D_xzz + alpha_ * Q_zz;
         F.D_yxx = -beta_X * D_yxx;
         F.D_yxy = -beta_X * D_yxy;
         F.D_yxz = -beta_X * D_yxz;
         F.D_yyy = -beta_X * D_yyy;
         F.D_yyz = -beta_X * D_yyz;
         F.D_yzz = -beta_X * D_yzz;
         F.D_zxx = -beta_X * D_zxx;
         F.D_zxy = -beta_X * D_zxy;
         F.D_zxz = -beta_X * D_zxz;
         F.D_zyy = -beta_X * D_zyy;
         F.D_zyz = -beta_X * D_zyz;
         F.D_zzz = -beta_X * D_zzz;
         
         F.Theta_ = -beta_X * Theta_ + alpha_ * (gamma_XX * V_x + gamma_XY * V_y + gamma_XZ * V_z);
         
         F.Z_x = -beta_X * Z_x - alpha_ * (K_Xx - K_ + Theta_);
         F.Z_y = -beta_X * Z_y - alpha_ * K_Xy;
         F.Z_z = -beta_X * Z_z - alpha_ * K_Xz;
         
         break;
      }
      
      case 1: {
         
         T const  D_yXX = gamma_XX * (D_yxx * gamma_XX + D_yxy * gamma_XY + D_yxz * gamma_XZ) + gamma_XY * (D_yxy * gamma_XX + D_yyy * gamma_XY + D_yyz * gamma_XZ) + gamma_XZ * (D_yxz * gamma_XX + D_yyz * gamma_XY + D_yzz * gamma_XZ);
         T const  D_yXY = gamma_XX * (D_yxx * gamma_XY + D_yxy * gamma_YY + D_yxz * gamma_YZ) + gamma_XY * (D_yxy * gamma_XY + D_yyy * gamma_YY + D_yyz * gamma_YZ) + gamma_XZ * (D_yxz * gamma_XY + D_yyz * gamma_YY + D_yzz * gamma_YZ);
         T const  D_yXZ = gamma_XX * (D_yxx * gamma_XZ + D_yxy * gamma_YZ + D_yxz * gamma_ZZ) + gamma_XY * (D_yxy * gamma_XZ + D_yyy * gamma_YZ + D_yyz * gamma_ZZ) + gamma_XZ * (D_yxz * gamma_XZ + D_yyz * gamma_YZ + D_yzz * gamma_ZZ);
         T const  D_yYY = gamma_XY * (D_yxx * gamma_XY + D_yxy * gamma_YY + D_yxz * gamma_YZ) + gamma_YY * (D_yxy * gamma_XY + D_yyy * gamma_YY + D_yyz * gamma_YZ) + gamma_YZ * (D_yxz * gamma_XY + D_yyz * gamma_YY + D_yzz * gamma_YZ);
         T const  D_yYZ = gamma_XY * (D_yxx * gamma_XZ + D_yxy * gamma_YZ + D_yxz * gamma_ZZ) + gamma_YY * (D_yxy * gamma_XZ + D_yyy * gamma_YZ + D_yyz * gamma_ZZ) + gamma_YZ * (D_yxz * gamma_XZ + D_yyz * gamma_YZ + D_yzz * gamma_ZZ);
         T const  D_yZZ = gamma_XZ * (D_yxx * gamma_XZ + D_yxy * gamma_YZ + D_yxz * gamma_ZZ) + gamma_YZ * (D_yxy * gamma_XZ + D_yyy * gamma_YZ + D_yyz * gamma_ZZ) + gamma_ZZ * (D_yxz * gamma_XZ + D_yyz * gamma_YZ + D_yzz * gamma_ZZ);
         
         T const  lambda_xxx = D_xxx + A_x + D_x - T(2) * (E_x + Z_x);
         T const  lambda_xxy = D_xxy + T(0.5) * (A_y + D_y - T(2) * (E_y + Z_y));
         T const  lambda_xxz = D_xxz + T(0.5) * (A_z + D_z - T(2) * (E_z + Z_z));
         T const &lambda_xyy = D_xyy;
         T const &lambda_xyz = D_xyz;
         T const &lambda_xzz = D_xzz;
         T const &lambda_yxx = D_yxx;
         T const  lambda_yxy = D_yxy + T(0.5) * (A_x + D_x - T(2) * (E_x + Z_x));
         T const &lambda_yxz = D_yxz;
         T const  lambda_yyy = D_yyy + A_y + D_y - T(2) * (E_y + Z_y);
         T const  lambda_yyz = D_yyz + T(0.5) * (A_z + D_z - T(2) * (E_z + Z_z));
         T const &lambda_yzz = D_yzz;
         T const &lambda_zxx = D_zxx;
         T const &lambda_zxy = D_zxy;
         T const  lambda_zxz = D_zxz + T(0.5) * (A_x + D_x - T(2) * (E_x + Z_x));
         T const &lambda_zyy = D_zyy;
         T const  lambda_zyz = D_zyz + T(0.5) * (A_y + D_y - T(2) * (E_y + Z_y));
         T const  lambda_zzz = D_zzz + A_z + D_z - T(2) * (E_z + Z_z);
         T const  lambda_Yxx = gamma_XY * lambda_xxx + gamma_YY * lambda_yxx + gamma_YZ * lambda_zxx;
         T const  lambda_Yxy = gamma_XY * lambda_xxy + gamma_YY * lambda_yxy + gamma_YZ * lambda_zxy;
         T const  lambda_Yxz = gamma_XY * lambda_xxz + gamma_YY * lambda_yxz + gamma_YZ * lambda_zxz;
         T const  lambda_Yyy = gamma_XY * lambda_xyy + gamma_YY * lambda_yyy + gamma_YZ * lambda_zyy;
         T const  lambda_Yyz = gamma_XY * lambda_xyz + gamma_YY * lambda_yyz + gamma_YZ * lambda_zyz;
         T const  lambda_Yzz = gamma_XY * lambda_xzz + gamma_YY * lambda_yzz + gamma_YZ * lambda_zzz;
         
         F.alpha_ = F.beta_X = F.beta_Y = F.beta_Z = F.gamma_xx = F.gamma_xy = F.gamma_xz = F.gamma_yy = F.gamma_yz = F.gamma_zz = T(0);
         
         F.K_xx = -beta_Y * K_xx + alpha_ * lambda_Yxx;
         F.K_xy = -beta_Y * K_xy + alpha_ * lambda_Yxy;
         F.K_xz = -beta_Y * K_xz + alpha_ * lambda_Yxz;
         F.K_yy = -beta_Y * K_yy + alpha_ * lambda_Yyy;
         F.K_yz = -beta_Y * K_yz + alpha_ * lambda_Yyz;
         F.K_zz = -beta_Y * K_zz + alpha_ * lambda_Yzz;
         
         F.A_x = -beta_Y * A_x;
         F.A_y = -beta_Y * A_y + alpha_ * Q_;
         F.A_z = -beta_Y * A_z;
         
         F.B_xX = -beta_Y * B_xX;
         F.B_xY = -beta_Y * B_xY;
         F.B_xZ = -beta_Y * B_xZ;
         F.B_yX = -beta_Y * B_yX + alpha_ * (gamma_XX * Q_x + gamma_XY * Q_y + gamma_XZ * Q_z);
         F.B_yY = -beta_Y * B_yY + alpha_ * (gamma_XY * Q_x + gamma_YY * Q_y + gamma_YZ * Q_z);
         F.B_yZ = -beta_Y * B_yZ + alpha_ * (gamma_XZ * Q_x + gamma_YZ * Q_y + gamma_ZZ * Q_z);
         F.B_zX = -beta_Y * B_zX;
         F.B_zY = -beta_Y * B_zY;
         F.B_zZ = -beta_Y * B_zZ;
         
         F.D_xxx = -beta_Y * D_xxx;
         F.D_xxy = -beta_Y * D_xxy;
         F.D_xxz = -beta_Y * D_xxz;
         F.D_xyy = -beta_Y * D_xyy;
         F.D_xyz = -beta_Y * D_xyz;
         F.D_xzz = -beta_Y * D_xzz;
         F.D_yxx = -beta_Y * D_yxx + alpha_ * Q_xx;
         F.D_yxy = -beta_Y * D_yxy + alpha_ * Q_xy;
         F.D_yxz = -beta_Y * D_yxz + alpha_ * Q_xz;
         F.D_yyy = -beta_Y * D_yyy + alpha_ * Q_yy;
         F.D_yyz = -beta_Y * D_yyz + alpha_ * Q_yz;
         F.D_yzz = -beta_Y * D_yzz + alpha_ * Q_zz;
         F.D_zxx = -beta_Y * D_zxx;
         F.D_zxy = -beta_Y * D_zxy;
         F.D_zxz = -beta_Y * D_zxz;
         F.D_zyy = -beta_Y * D_zyy;
         F.D_zyz = -beta_Y * D_zyz;
         F.D_zzz = -beta_Y * D_zzz;
         
         F.Theta_ = -beta_Y * Theta_ + alpha_ * (gamma_XY * V_x + gamma_YY * V_y + gamma_YZ * V_z);
         
         F.Z_x = -beta_Y * Z_x - alpha_ * K_Yx;
         F.Z_y = -beta_Y * Z_y - alpha_ * (K_Yy - K_ + Theta_);
         F.Z_z = -beta_Y * Z_z - alpha_ * K_Yz;
         
         break;
      }
      
      case 2: {
         
         T const  D_zXX = gamma_XX * (D_zxx * gamma_XX + D_zxy * gamma_XY + D_zxz * gamma_XZ) + gamma_XY * (D_zxy * gamma_XX + D_zyy * gamma_XY + D_zyz * gamma_XZ) + gamma_XZ * (D_zxz * gamma_XX + D_zyz * gamma_XY + D_zzz * gamma_XZ);
         T const  D_zXY = gamma_XX * (D_zxx * gamma_XY + D_zxy * gamma_YY + D_zxz * gamma_YZ) + gamma_XY * (D_zxy * gamma_XY + D_zyy * gamma_YY + D_zyz * gamma_YZ) + gamma_XZ * (D_zxz * gamma_XY + D_zyz * gamma_YY + D_zzz * gamma_YZ);
         T const  D_zXZ = gamma_XX * (D_zxx * gamma_XZ + D_zxy * gamma_YZ + D_zxz * gamma_ZZ) + gamma_XY * (D_zxy * gamma_XZ + D_zyy * gamma_YZ + D_zyz * gamma_ZZ) + gamma_XZ * (D_zxz * gamma_XZ + D_zyz * gamma_YZ + D_zzz * gamma_ZZ);
         T const  D_zYY = gamma_XY * (D_zxx * gamma_XY + D_zxy * gamma_YY + D_zxz * gamma_YZ) + gamma_YY * (D_zxy * gamma_XY + D_zyy * gamma_YY + D_zyz * gamma_YZ) + gamma_YZ * (D_zxz * gamma_XY + D_zyz * gamma_YY + D_zzz * gamma_YZ);
         T const  D_zYZ = gamma_XY * (D_zxx * gamma_XZ + D_zxy * gamma_YZ + D_zxz * gamma_ZZ) + gamma_YY * (D_zxy * gamma_XZ + D_zyy * gamma_YZ + D_zyz * gamma_ZZ) + gamma_YZ * (D_zxz * gamma_XZ + D_zyz * gamma_YZ + D_zzz * gamma_ZZ);
         T const  D_zZZ = gamma_XZ * (D_zxx * gamma_XZ + D_zxy * gamma_YZ + D_zxz * gamma_ZZ) + gamma_YZ * (D_zxy * gamma_XZ + D_zyy * gamma_YZ + D_zyz * gamma_ZZ) + gamma_ZZ * (D_zxz * gamma_XZ + D_zyz * gamma_YZ + D_zzz * gamma_ZZ);
         
         T const  lambda_xxx = D_xxx + A_x + D_x - T(2) * (E_x + Z_x);
         T const  lambda_xxy = D_xxy + T(0.5) * (A_y + D_y - T(2) * (E_y + Z_y));
         T const  lambda_xxz = D_xxz + T(0.5) * (A_z + D_z - T(2) * (E_z + Z_z));
         T const &lambda_xyy = D_xyy;
         T const &lambda_xyz = D_xyz;
         T const &lambda_xzz = D_xzz;
         T const &lambda_yxx = D_yxx;
         T const  lambda_yxy = D_yxy + T(0.5) * (A_x + D_x - T(2) * (E_x + Z_x));
         T const &lambda_yxz = D_yxz;
         T const  lambda_yyy = D_yyy + A_y + D_y - T(2) * (E_y + Z_y);
         T const  lambda_yyz = D_yyz + T(0.5) * (A_z + D_z - T(2) * (E_z + Z_z));
         T const &lambda_yzz = D_yzz;
         T const &lambda_zxx = D_zxx;
         T const &lambda_zxy = D_zxy;
         T const  lambda_zxz = D_zxz + T(0.5) * (A_x + D_x - T(2) * (E_x + Z_x));
         T const &lambda_zyy = D_zyy;
         T const  lambda_zyz = D_zyz + T(0.5) * (A_y + D_y - T(2) * (E_y + Z_y));
         T const  lambda_zzz = D_zzz + A_z + D_z - T(2) * (E_z + Z_z);
         T const  lambda_Zxx = gamma_XZ * lambda_xxx + gamma_YZ * lambda_yxx + gamma_ZZ * lambda_zxx;
         T const  lambda_Zxy = gamma_XZ * lambda_xxy + gamma_YZ * lambda_yxy + gamma_ZZ * lambda_zxy;
         T const  lambda_Zxz = gamma_XZ * lambda_xxz + gamma_YZ * lambda_yxz + gamma_ZZ * lambda_zxz;
         T const  lambda_Zyy = gamma_XZ * lambda_xyy + gamma_YZ * lambda_yyy + gamma_ZZ * lambda_zyy;
         T const  lambda_Zyz = gamma_XZ * lambda_xyz + gamma_YZ * lambda_yyz + gamma_ZZ * lambda_zyz;
         T const  lambda_Zzz = gamma_XZ * lambda_xzz + gamma_YZ * lambda_yzz + gamma_ZZ * lambda_zzz;
         
         F.alpha_ = F.beta_X = F.beta_Y = F.beta_Z = F.gamma_xx = F.gamma_xy = F.gamma_xz = F.gamma_yy = F.gamma_yz = F.gamma_zz = T(0);
         
         F.K_xx = -beta_Z * K_xx + alpha_ * lambda_Zxx;
         F.K_xy = -beta_Z * K_xy + alpha_ * lambda_Zxy;
         F.K_xz = -beta_Z * K_xz + alpha_ * lambda_Zxz;
         F.K_yy = -beta_Z * K_yy + alpha_ * lambda_Zyy;
         F.K_yz = -beta_Z * K_yz + alpha_ * lambda_Zyz;
         F.K_zz = -beta_Z * K_zz + alpha_ * lambda_Zzz;
         
         F.A_x = -beta_Z * A_x;
         F.A_y = -beta_Z * A_y;
         F.A_z = -beta_Z * A_z + alpha_ * Q_;
         
         F.B_xX = -beta_Z * B_xX;
         F.B_xY = -beta_Z * B_xY;
         F.B_xZ = -beta_Z * B_xZ;
         F.B_yX = -beta_Z * B_yX;
         F.B_yY = -beta_Z * B_yY;
         F.B_yZ = -beta_Z * B_yZ;
         F.B_zX = -beta_Z * B_zX + alpha_ * (gamma_XX * Q_x + gamma_XY * Q_y + gamma_XZ * Q_z);
         F.B_zY = -beta_Z * B_zY + alpha_ * (gamma_XY * Q_x + gamma_YY * Q_y + gamma_YZ * Q_z);
         F.B_zZ = -beta_Z * B_zZ + alpha_ * (gamma_XZ * Q_x + gamma_YZ * Q_y + gamma_ZZ * Q_z);
         
         F.D_xxx = -beta_Z * D_xxx;
         F.D_xxy = -beta_Z * D_xxy;
         F.D_xxz = -beta_Z * D_xxz;
         F.D_xyy = -beta_Z * D_xyy;
         F.D_xyz = -beta_Z * D_xyz;
         F.D_xzz = -beta_Z * D_xzz;
         F.D_yxx = -beta_Z * D_yxx;
         F.D_yxy = -beta_Z * D_yxy;
         F.D_yxz = -beta_Z * D_yxz;
         F.D_yyy = -beta_Z * D_yyy;
         F.D_yyz = -beta_Z * D_yyz;
         F.D_yzz = -beta_Z * D_yzz;
         F.D_zxx = -beta_Z * D_zxx + alpha_ * Q_xx;
         F.D_zxy = -beta_Z * D_zxy + alpha_ * Q_xy;
         F.D_zxz = -beta_Z * D_zxz + alpha_ * Q_xz;
         F.D_zyy = -beta_Z * D_zyy + alpha_ * Q_yy;
         F.D_zyz = -beta_Z * D_zyz + alpha_ * Q_yz;
         F.D_zzz = -beta_Z * D_zzz + alpha_ * Q_zz;
         
         F.Theta_ = -beta_Z * Theta_ + alpha_ * (gamma_XZ * V_x + gamma_YZ * V_y + gamma_ZZ * V_z);
         
         F.Z_x = -beta_Z * Z_x - alpha_ * K_Zx;
         F.Z_y = -beta_Z * Z_y - alpha_ * K_Zy;
         F.Z_z = -beta_Z * Z_z - alpha_ * (K_Zz - K_ + Theta_);
         
         break;
      }
      
      default:
         throw std::invalid_argument( "Index for Spacetime flux must be 0, 1, or 2." );
      
      }
      
      return F;
   }
   
   Spacetime source() const {
      Spacetime S;
      
      T gamma_XX = gamma_yy * gamma_zz - gamma_yz * gamma_yz;
      T gamma_YY = gamma_xx * gamma_zz - gamma_xz * gamma_xz;
      T gamma_ZZ = gamma_xx * gamma_yy - gamma_xy * gamma_xy;
      T gamma_XY = gamma_zz * gamma_xy - gamma_xz * gamma_yz;
      T gamma_XZ = gamma_xy * gamma_yz - gamma_yy * gamma_xz;
      T gamma_YZ = gamma_yz * gamma_xx - gamma_xz * gamma_xy;
      
      T const  gamma_ = gamma_xx * gamma_XX + gamma_xy * gamma_XY + gamma_xz * gamma_XZ;
      
      T const  igamma_ = T(1) / gamma_;
      
      gamma_XX *= igamma_;
      gamma_YY *= igamma_;
      gamma_ZZ *= igamma_;
      gamma_XY *= igamma_;
      gamma_XZ *= igamma_;
      gamma_YZ *= igamma_;
      
      T const  B_ = B_xX + B_yY + B_zZ;
      
      T const  K_ = gamma_XX * K_xx + gamma_YY * K_yy + gamma_ZZ * K_zz + T(2) * (gamma_XY * K_xy + gamma_XZ * K_xz + gamma_YZ * K_yz);
      
      T const  D_x = D_xxx * gamma_XX + D_xyy * gamma_YY + D_xzz * gamma_ZZ + T(2) * (D_xxy * gamma_XY + D_xxz * gamma_XZ + D_xyz * gamma_YZ);
      T const  D_y = D_yxx * gamma_XX + D_yyy * gamma_YY + D_yzz * gamma_ZZ + T(2) * (D_yxy * gamma_XY + D_yxz * gamma_XZ + D_yyz * gamma_YZ);
      T const  D_z = D_zxx * gamma_XX + D_zyy * gamma_YY + D_zzz * gamma_ZZ + T(2) * (D_zxy * gamma_XY + D_zxz * gamma_XZ + D_zyz * gamma_YZ);
      
      T const  E_x = D_xxx * gamma_XX + D_xxy * gamma_XY + D_xxz * gamma_XZ + D_yxx * gamma_XY + D_yxy * gamma_YY + D_yxz * gamma_YZ + D_zxx * gamma_XZ + D_zxy * gamma_YZ + D_zxz * gamma_ZZ;
      T const  E_y = D_xxy * gamma_XX + D_xyy * gamma_XY + D_xyz * gamma_XZ + D_yxy * gamma_XY + D_yyy * gamma_YY + D_yyz * gamma_YZ + D_zxy * gamma_XZ + D_zyy * gamma_YZ + D_zyz * gamma_ZZ;
      T const  E_z = D_xxz * gamma_XX + D_xyz * gamma_XY + D_xzz * gamma_XZ + D_yxz * gamma_XY + D_yyz * gamma_YY + D_yzz * gamma_YZ + D_zxz * gamma_XZ + D_zyz * gamma_YZ + D_zzz * gamma_ZZ;
      
      T const  D_X = D_x * gamma_XX + D_y * gamma_XY + D_z * gamma_XZ;
      T const  D_Y = D_x * gamma_XY + D_y * gamma_YY + D_z * gamma_YZ;
      T const  D_Z = D_x * gamma_XZ + D_y * gamma_YZ + D_z * gamma_ZZ;
      
      T const  E_X = E_x * gamma_XX + E_y * gamma_XY + E_z * gamma_XZ;
      T const  E_Y = E_x * gamma_XY + E_y * gamma_YY + E_z * gamma_YZ;
      T const  E_Z = E_x * gamma_XZ + E_y * gamma_YZ + E_z * gamma_ZZ;
      
      T const  Z_X = Z_x * gamma_XX + Z_y * gamma_XY + Z_z * gamma_XZ;
      T const  Z_Y = Z_x * gamma_XY + Z_y * gamma_YY + Z_z * gamma_YZ;
      T const  Z_Z = Z_x * gamma_XZ + Z_y * gamma_YZ + Z_z * gamma_ZZ;
      
      T const  K_Xx = gamma_XX * K_xx + gamma_XY * K_xy + gamma_XZ * K_xz;
      T const  K_Xy = gamma_XX * K_xy + gamma_XY * K_yy + gamma_XZ * K_yz;
      T const  K_Xz = gamma_XX * K_xz + gamma_XY * K_yz + gamma_XZ * K_zz;
      T const  K_Yx = gamma_XY * K_xx + gamma_YY * K_xy + gamma_YZ * K_xz;
      T const  K_Yy = gamma_XY * K_xy + gamma_YY * K_yy + gamma_YZ * K_yz;
      T const  K_Yz = gamma_XY * K_xz + gamma_YY * K_yz + gamma_YZ * K_zz;
      T const  K_Zx = gamma_XZ * K_xx + gamma_YZ * K_xy + gamma_ZZ * K_xz;
      T const  K_Zy = gamma_XZ * K_xy + gamma_YZ * K_yy + gamma_ZZ * K_yz;
      T const  K_Zz = gamma_XZ * K_xz + gamma_YZ * K_yz + gamma_ZZ * K_zz;
      
      T const  D_xXX = gamma_XX * (D_xxx * gamma_XX + D_xxy * gamma_XY + D_xxz * gamma_XZ) + gamma_XY * (D_xxy * gamma_XX + D_xyy * gamma_XY + D_xyz * gamma_XZ) + gamma_XZ * (D_xxz * gamma_XX + D_xyz * gamma_XY + D_xzz * gamma_XZ);
      T const  D_xXY = gamma_XX * (D_xxx * gamma_XY + D_xxy * gamma_YY + D_xxz * gamma_YZ) + gamma_XY * (D_xxy * gamma_XY + D_xyy * gamma_YY + D_xyz * gamma_YZ) + gamma_XZ * (D_xxz * gamma_XY + D_xyz * gamma_YY + D_xzz * gamma_YZ);
      T const  D_xXZ = gamma_XX * (D_xxx * gamma_XZ + D_xxy * gamma_YZ + D_xxz * gamma_ZZ) + gamma_XY * (D_xxy * gamma_XZ + D_xyy * gamma_YZ + D_xyz * gamma_ZZ) + gamma_XZ * (D_xxz * gamma_XZ + D_xyz * gamma_YZ + D_xzz * gamma_ZZ);
      T const  D_xYY = gamma_XY * (D_xxx * gamma_XY + D_xxy * gamma_YY + D_xxz * gamma_YZ) + gamma_YY * (D_xxy * gamma_XY + D_xyy * gamma_YY + D_xyz * gamma_YZ) + gamma_YZ * (D_xxz * gamma_XY + D_xyz * gamma_YY + D_xzz * gamma_YZ);
      T const  D_xYZ = gamma_XY * (D_xxx * gamma_XZ + D_xxy * gamma_YZ + D_xxz * gamma_ZZ) + gamma_YY * (D_xxy * gamma_XZ + D_xyy * gamma_YZ + D_xyz * gamma_ZZ) + gamma_YZ * (D_xxz * gamma_XZ + D_xyz * gamma_YZ + D_xzz * gamma_ZZ);
      T const  D_xZZ = gamma_XZ * (D_xxx * gamma_XZ + D_xxy * gamma_YZ + D_xxz * gamma_ZZ) + gamma_YZ * (D_xxy * gamma_XZ + D_xyy * gamma_YZ + D_xyz * gamma_ZZ) + gamma_ZZ * (D_xxz * gamma_XZ + D_xyz * gamma_YZ + D_xzz * gamma_ZZ);
      T const  D_yXX = gamma_XX * (D_yxx * gamma_XX + D_yxy * gamma_XY + D_yxz * gamma_XZ) + gamma_XY * (D_yxy * gamma_XX + D_yyy * gamma_XY + D_yyz * gamma_XZ) + gamma_XZ * (D_yxz * gamma_XX + D_yyz * gamma_XY + D_yzz * gamma_XZ);
      T const  D_yXY = gamma_XX * (D_yxx * gamma_XY + D_yxy * gamma_YY + D_yxz * gamma_YZ) + gamma_XY * (D_yxy * gamma_XY + D_yyy * gamma_YY + D_yyz * gamma_YZ) + gamma_XZ * (D_yxz * gamma_XY + D_yyz * gamma_YY + D_yzz * gamma_YZ);
      T const  D_yXZ = gamma_XX * (D_yxx * gamma_XZ + D_yxy * gamma_YZ + D_yxz * gamma_ZZ) + gamma_XY * (D_yxy * gamma_XZ + D_yyy * gamma_YZ + D_yyz * gamma_ZZ) + gamma_XZ * (D_yxz * gamma_XZ + D_yyz * gamma_YZ + D_yzz * gamma_ZZ);
      T const  D_yYY = gamma_XY * (D_yxx * gamma_XY + D_yxy * gamma_YY + D_yxz * gamma_YZ) + gamma_YY * (D_yxy * gamma_XY + D_yyy * gamma_YY + D_yyz * gamma_YZ) + gamma_YZ * (D_yxz * gamma_XY + D_yyz * gamma_YY + D_yzz * gamma_YZ);
      T const  D_yYZ = gamma_XY * (D_yxx * gamma_XZ + D_yxy * gamma_YZ + D_yxz * gamma_ZZ) + gamma_YY * (D_yxy * gamma_XZ + D_yyy * gamma_YZ + D_yyz * gamma_ZZ) + gamma_YZ * (D_yxz * gamma_XZ + D_yyz * gamma_YZ + D_yzz * gamma_ZZ);
      T const  D_yZZ = gamma_XZ * (D_yxx * gamma_XZ + D_yxy * gamma_YZ + D_yxz * gamma_ZZ) + gamma_YZ * (D_yxy * gamma_XZ + D_yyy * gamma_YZ + D_yyz * gamma_ZZ) + gamma_ZZ * (D_yxz * gamma_XZ + D_yyz * gamma_YZ + D_yzz * gamma_ZZ);
      T const  D_zXX = gamma_XX * (D_zxx * gamma_XX + D_zxy * gamma_XY + D_zxz * gamma_XZ) + gamma_XY * (D_zxy * gamma_XX + D_zyy * gamma_XY + D_zyz * gamma_XZ) + gamma_XZ * (D_zxz * gamma_XX + D_zyz * gamma_XY + D_zzz * gamma_XZ);
      T const  D_zXY = gamma_XX * (D_zxx * gamma_XY + D_zxy * gamma_YY + D_zxz * gamma_YZ) + gamma_XY * (D_zxy * gamma_XY + D_zyy * gamma_YY + D_zyz * gamma_YZ) + gamma_XZ * (D_zxz * gamma_XY + D_zyz * gamma_YY + D_zzz * gamma_YZ);
      T const  D_zXZ = gamma_XX * (D_zxx * gamma_XZ + D_zxy * gamma_YZ + D_zxz * gamma_ZZ) + gamma_XY * (D_zxy * gamma_XZ + D_zyy * gamma_YZ + D_zyz * gamma_ZZ) + gamma_XZ * (D_zxz * gamma_XZ + D_zyz * gamma_YZ + D_zzz * gamma_ZZ);
      T const  D_zYY = gamma_XY * (D_zxx * gamma_XY + D_zxy * gamma_YY + D_zxz * gamma_YZ) + gamma_YY * (D_zxy * gamma_XY + D_zyy * gamma_YY + D_zyz * gamma_YZ) + gamma_YZ * (D_zxz * gamma_XY + D_zyz * gamma_YY + D_zzz * gamma_YZ);
      T const  D_zYZ = gamma_XY * (D_zxx * gamma_XZ + D_zxy * gamma_YZ + D_zxz * gamma_ZZ) + gamma_YY * (D_zxy * gamma_XZ + D_zyy * gamma_YZ + D_zyz * gamma_ZZ) + gamma_YZ * (D_zxz * gamma_XZ + D_zyz * gamma_YZ + D_zzz * gamma_ZZ);
      T const  D_zZZ = gamma_XZ * (D_zxx * gamma_XZ + D_zxy * gamma_YZ + D_zxz * gamma_ZZ) + gamma_YZ * (D_zxy * gamma_XZ + D_zyy * gamma_YZ + D_zyz * gamma_ZZ) + gamma_ZZ * (D_zxz * gamma_XZ + D_zyz * gamma_YZ + D_zzz * gamma_ZZ);
      
      T const &Gamma_xxx = D_xxx;
      T const  Gamma_xyy = T(2) * D_yxy - D_xyy;
      T const  Gamma_xzz = T(2) * D_zxz - D_xzz;
      T const &Gamma_xxy = D_yxx;
      T const &Gamma_xxz = D_zxx;
      T const  Gamma_xyz = D_yxz + D_zxy - D_xyz;
      T const  Gamma_yxx = T(2) * D_xxy - D_yxx;
      T const &Gamma_yyy = D_yyy;
      T const  Gamma_yzz = T(2) * D_zyz - D_yzz;
      T const &Gamma_yxy = D_yxy;
      T const  Gamma_yxz = D_xyz + D_zxy - D_yxz;
      T const &Gamma_yyz = D_zyy;
      T const  Gamma_zxx = T(2) * D_xxz - D_zxx;
      T const  Gamma_zyy = T(2) * D_yyz - D_zyy;
      T const &Gamma_zzz = D_zzz;
      T const  Gamma_zxy = D_xyz + D_yxz - D_zxy;
      T const &Gamma_zxz = D_zxz;
      T const &Gamma_zyz = D_zyz;
      
      T const  Gamma_Xxx = T(0.5) * (gamma_XX * Gamma_xxx + gamma_XY * Gamma_yxx + gamma_XZ * Gamma_zxx);
      T const  Gamma_Xyy = T(0.5) * (gamma_XX * Gamma_xyy + gamma_XY * Gamma_yyy + gamma_XZ * Gamma_zyy);
      T const  Gamma_Xzz = T(0.5) * (gamma_XX * Gamma_xzz + gamma_XY * Gamma_yzz + gamma_XZ * Gamma_zzz);
      T const  Gamma_Xxy = T(0.5) * (gamma_XX * Gamma_xxy + gamma_XY * Gamma_yxy + gamma_XZ * Gamma_zxy);
      T const  Gamma_Xxz = T(0.5) * (gamma_XX * Gamma_xxz + gamma_XY * Gamma_yxz + gamma_XZ * Gamma_zxz);
      T const  Gamma_Xyz = T(0.5) * (gamma_XX * Gamma_xyz + gamma_XY * Gamma_yyz + gamma_XZ * Gamma_zyz);
      T const  Gamma_Yxx = T(0.5) * (gamma_XY * Gamma_xxx + gamma_YY * Gamma_yxx + gamma_YZ * Gamma_zxx);
      T const  Gamma_Yyy = T(0.5) * (gamma_XY * Gamma_xyy + gamma_YY * Gamma_yyy + gamma_YZ * Gamma_zyy);
      T const  Gamma_Yzz = T(0.5) * (gamma_XY * Gamma_xzz + gamma_YY * Gamma_yzz + gamma_YZ * Gamma_zzz);
      T const  Gamma_Yxy = T(0.5) * (gamma_XY * Gamma_xxy + gamma_YY * Gamma_yxy + gamma_YZ * Gamma_zxy);
      T const  Gamma_Yxz = T(0.5) * (gamma_XY * Gamma_xxz + gamma_YY * Gamma_yxz + gamma_YZ * Gamma_zxz);
      T const  Gamma_Yyz = T(0.5) * (gamma_XY * Gamma_xyz + gamma_YY * Gamma_yyz + gamma_YZ * Gamma_zyz);
      T const  Gamma_Zxx = T(0.5) * (gamma_XZ * Gamma_xxx + gamma_YZ * Gamma_yxx + gamma_ZZ * Gamma_zxx);
      T const  Gamma_Zyy = T(0.5) * (gamma_XZ * Gamma_xyy + gamma_YZ * Gamma_yyy + gamma_ZZ * Gamma_zyy);
      T const  Gamma_Zzz = T(0.5) * (gamma_XZ * Gamma_xzz + gamma_YZ * Gamma_yzz + gamma_ZZ * Gamma_zzz);
      T const  Gamma_Zxy = T(0.5) * (gamma_XZ * Gamma_xxy + gamma_YZ * Gamma_yxy + gamma_ZZ * Gamma_zxy);
      T const  Gamma_Zxz = T(0.5) * (gamma_XZ * Gamma_xxz + gamma_YZ * Gamma_yxz + gamma_ZZ * Gamma_zxz);
      T const  Gamma_Zyz = T(0.5) * (gamma_XZ * Gamma_xyz + gamma_YZ * Gamma_yyz + gamma_ZZ * Gamma_zyz);
      
      T const  V_x = D_x - E_x - Z_x;
      T const  V_y = D_y - E_y - Z_y;
      T const  V_z = D_z - E_z - Z_z;
      
      T const  ahpla_ = T(1) / alpha_;
      
      T const  Q_ = K_ - T(2) * Theta_;
      
      T const  Q_x = alpha_ * (A_x - D_x + T(2) * V_x);
      T const  Q_y = alpha_ * (A_y - D_y + T(2) * V_y);
      T const  Q_z = alpha_ * (A_z - D_z + T(2) * V_z);
      
      T const  Q_xx = K_xx - ahpla_ * B_xX;
      T const  Q_yy = K_yy - ahpla_ * B_yY;
      T const  Q_zz = K_zz - ahpla_ * B_zZ;
      T const  Q_xy = K_xy - T(0.5) * ahpla_ * (B_xY + B_yX);
      T const  Q_xz = K_xz - T(0.5) * ahpla_ * (B_xZ + B_zX);
      T const  Q_yz = K_yz - T(0.5) * ahpla_ * (B_yZ + B_zY);
      
      S.alpha_ = alpha_ * (beta_X * A_x + beta_Y * A_y + beta_Z * A_z - alpha_ * Q_);
      
      S.beta_X = beta_X * B_xX + beta_Y * B_yX + beta_Z * B_zX - alpha_ * (gamma_XX * Q_x + gamma_XY * Q_y + gamma_XZ * Q_z);
      S.beta_Y = beta_X * B_xY + beta_Y * B_yY + beta_Z * B_zY - alpha_ * (gamma_XY * Q_x + gamma_YY * Q_y + gamma_YZ * Q_z);
      S.beta_Z = beta_X * B_xZ + beta_Y * B_yZ + beta_Z * B_zZ - alpha_ * (gamma_XZ * Q_x + gamma_YZ * Q_y + gamma_ZZ * Q_z);
      
      S.gamma_xx = T(2) * (beta_X * D_xxx + beta_Y * D_yxx + beta_Z * D_zxx - alpha_ * K_xx + B_xX);
      S.gamma_yy = T(2) * (beta_X * D_xyy + beta_Y * D_yyy + beta_Z * D_zyy - alpha_ * K_yy + B_yY);
      S.gamma_zz = T(2) * (beta_X * D_xzz + beta_Y * D_yzz + beta_Z * D_zzz - alpha_ * K_zz + B_zZ);
      S.gamma_xy = T(2) * (beta_X * D_xxy + beta_Y * D_yxy + beta_Z * D_zxy - alpha_ * K_xy) + B_xY + B_yX;
      S.gamma_xz = T(2) * (beta_X * D_xxz + beta_Y * D_yxz + beta_Z * D_zxz - alpha_ * K_xz) + B_xZ + B_zX;
      S.gamma_yz = T(2) * (beta_X * D_xyz + beta_Y * D_yyz + beta_Z * D_zyz - alpha_ * K_yz) + B_yZ + B_zY;
      
      S.K_xx = -B_ * K_xx + K_xx * B_xX + K_xx * B_xX + K_xy * B_xY + K_xy * B_xY + K_xz * B_xZ + K_xz * B_xZ +  alpha_ * (
         (A_x - T(2) * Z_x) * D_x + (D_x - T(2) * Z_x) * Gamma_Xxx + (D_y - T(2) * Z_y) * Gamma_Yxx + (D_z - T(2) * Z_z) * Gamma_Zxx +
         -(Gamma_Xxx * Gamma_Xxx + Gamma_Yxx * Gamma_Xxy + Gamma_Zxx * Gamma_Xxz + Gamma_Xxy * Gamma_Yxx + Gamma_Yxy * Gamma_Yxy + Gamma_Zxy * Gamma_Yxz + Gamma_Xxz * Gamma_Zxx + Gamma_Yxz * Gamma_Zxy + Gamma_Zxz * Gamma_Zxz) +
         -T(2) * (K_Xx * K_xx + K_Yx * K_xy + K_Zx * K_xz + K_Xx * K_xx + K_Yx * K_xy + K_Zx * K_xz + K_Xx * K_xx + K_Yx * K_xy + K_Zx * K_xz) +
         Q_ * K_xx);
      S.K_yy = -B_ * K_yy + K_xy * B_yX + K_xy * B_yX + K_yy * B_yY + K_yy * B_yY + K_yz * B_yZ + K_yz * B_yZ +  alpha_ * (
         (A_y - T(2) * Z_y) * D_y + (D_x - T(2) * Z_x) * Gamma_Xyy + (D_y - T(2) * Z_y) * Gamma_Yyy + (D_z - T(2) * Z_z) * Gamma_Zyy +
         -(Gamma_Xxy * Gamma_Xxy + Gamma_Yxy * Gamma_Xyy + Gamma_Zxy * Gamma_Xyz + Gamma_Xyy * Gamma_Yxy + Gamma_Yyy * Gamma_Yyy + Gamma_Zyy * Gamma_Yyz + Gamma_Xyz * Gamma_Zxy + Gamma_Yyz * Gamma_Zyy + Gamma_Zyz * Gamma_Zyz) +
         -T(2) * (K_Xy * K_xy + K_Yy * K_yy + K_Zy * K_yz + K_Xy * K_xy + K_Yy * K_yy + K_Zy * K_yz + K_Xy * K_xy + K_Yy * K_yy + K_Zy * K_yz) +
         Q_ * K_yy);
      S.K_zz = -B_ * K_zz + K_xz * B_zX + K_xz * B_zX + K_yz * B_zY + K_yz * B_zY + K_zz * B_zZ + K_zz * B_zZ +  alpha_ * (
         (A_z - T(2) * Z_z) * D_z + (D_x - T(2) * Z_x) * Gamma_Xzz + (D_y - T(2) * Z_y) * Gamma_Yzz + (D_z - T(2) * Z_z) * Gamma_Zzz +
         -(Gamma_Xxz * Gamma_Xxz + Gamma_Yxz * Gamma_Xyz + Gamma_Zxz * Gamma_Xzz + Gamma_Xyz * Gamma_Yxz + Gamma_Yyz * Gamma_Yyz + Gamma_Zyz * Gamma_Yzz + Gamma_Xzz * Gamma_Zxz + Gamma_Yzz * Gamma_Zyz + Gamma_Zzz * Gamma_Zzz) +
         -T(2) * (K_Xz * K_xz + K_Yz * K_yz + K_Zz * K_zz + K_Xz * K_xz + K_Yz * K_yz + K_Zz * K_zz + K_Xz * K_xz + K_Yz * K_yz + K_Zz * K_zz) +
         Q_ * K_zz);
      S.K_xy = -B_ * K_xy + K_xx * B_yX + K_xy * B_xX + K_xy * B_yY + K_yy * B_xY + K_xz * B_yZ + K_yz * B_xZ +  alpha_ * (
         T(0.5) * ((A_x - T(2) * Z_x) * D_y + (A_y - T(2) * Z_y) * D_x) + (D_x - T(2) * Z_x) * Gamma_Xxy + (D_y - T(2) * Z_y) * Gamma_Yxy + (D_z - T(2) * Z_z) * Gamma_Zxy +
         -(Gamma_Xxy * Gamma_Xxx + Gamma_Yxy * Gamma_Xxy + Gamma_Zxy * Gamma_Xxz + Gamma_Xyy * Gamma_Yxx + Gamma_Yyy * Gamma_Yxy + Gamma_Zyy * Gamma_Yxz + Gamma_Xyz * Gamma_Zxx + Gamma_Yyz * Gamma_Zxy + Gamma_Zyz * Gamma_Zxz) +
         -T(2) * (K_Xx * K_xy + K_Yx * K_yy + K_Zx * K_yz + K_Xx * K_xy + K_Yx * K_yy + K_Zx * K_yz + K_Xx * K_xy + K_Yx * K_yy + K_Zx * K_yz) +
         Q_ * K_xy);
      S.K_xz = -B_ * K_xz + K_xx * B_zX + K_xz * B_xX + K_xy * B_zY + K_yz * B_xY + K_xz * B_zZ + K_zz * B_xZ +  alpha_ * (
         T(0.5) * ((A_x - T(2) * Z_x) * D_z + (A_z - T(2) * Z_z) * D_x) + (D_x - T(2) * Z_x) * Gamma_Xxz + (D_y - T(2) * Z_y) * Gamma_Yxz + (D_z - T(2) * Z_z) * Gamma_Zxz +
         -(Gamma_Xxz * Gamma_Xxx + Gamma_Yxz * Gamma_Xxy + Gamma_Zxz * Gamma_Xxz + Gamma_Xyz * Gamma_Yxx + Gamma_Yyz * Gamma_Yxy + Gamma_Zyz * Gamma_Yxz + Gamma_Xzz * Gamma_Zxx + Gamma_Yzz * Gamma_Zxy + Gamma_Zzz * Gamma_Zxz) +
         -T(2) * (K_Xx * K_xz + K_Yx * K_yz + K_Zx * K_zz + K_Xx * K_xz + K_Yx * K_yz + K_Zx * K_zz + K_Xx * K_xz + K_Yx * K_yz + K_Zx * K_zz) +
         Q_ * K_xz);
      S.K_yz = -B_ * K_yz + K_xy * B_zX + K_xz * B_yX + K_yy * B_zY + K_yz * B_yY + K_yz * B_zZ + K_zz * B_yZ +  alpha_ * (
         T(0.5) * ((A_y - T(2) * Z_y) * D_z + (A_z - T(2) * Z_z) * D_y) + (D_x - T(2) * Z_x) * Gamma_Xyz + (D_y - T(2) * Z_y) * Gamma_Yyz + (D_z - T(2) * Z_z) * Gamma_Zyz +
         -(Gamma_Xxz * Gamma_Xxy + Gamma_Yxz * Gamma_Xyy + Gamma_Zxz * Gamma_Xyz + Gamma_Xyz * Gamma_Yxy + Gamma_Yyz * Gamma_Yyy + Gamma_Zyz * Gamma_Yyz + Gamma_Xzz * Gamma_Zxy + Gamma_Yzz * Gamma_Zyy + Gamma_Zzz * Gamma_Zyz) +
         -T(2) * (K_Xy * K_xz + K_Yy * K_yz + K_Zy * K_zz + K_Xy * K_xz + K_Yy * K_yz + K_Zy * K_zz + K_Xy * K_xz + K_Yy * K_yz + K_Zy * K_zz) +
         Q_ * K_yz);
      
      S.A_x = B_xX * A_x + B_xY * A_y + B_xZ * A_z - B_ * A_x;
      S.A_y = B_yX * A_x + B_yY * A_y + B_yZ * A_z - B_ * A_y;
      S.A_z = B_zX * A_x + B_zY * A_y + B_zZ * A_z - B_ * A_z;
      
      S.B_xX = B_xX * B_xX + B_xY * B_yX + B_xZ * B_zX - B_ * B_xX;
      S.B_xY = B_xX * B_xY + B_xY * B_yY + B_xZ * B_zY - B_ * B_xY;
      S.B_xZ = B_xX * B_xZ + B_xY * B_yZ + B_xZ * B_zZ - B_ * B_xZ;
      S.B_yX = B_yX * B_xX + B_yY * B_yX + B_yZ * B_zX - B_ * B_yX;
      S.B_yY = B_yX * B_xY + B_yY * B_yY + B_yZ * B_zY - B_ * B_yY;
      S.B_yZ = B_yX * B_xZ + B_yY * B_yZ + B_yZ * B_zZ - B_ * B_yZ;
      S.B_zX = B_zX * B_xX + B_zY * B_yX + B_zZ * B_zX - B_ * B_zX;
      S.B_zY = B_zX * B_xY + B_zY * B_yY + B_zZ * B_zY - B_ * B_zY;
      S.B_zZ = B_zX * B_xZ + B_zY * B_yZ + B_zZ * B_zZ - B_ * B_zZ;
      
      S.D_xxx = B_xX * D_xxx + B_xY * D_yxx + B_xZ * D_zxx - B_ * D_xxx;
      S.D_xxy = B_xX * D_xxy + B_xY * D_yxy + B_xZ * D_zxy - B_ * D_xxy;
      S.D_xxz = B_xX * D_xxz + B_xY * D_yxz + B_xZ * D_zxz - B_ * D_xxz;
      S.D_xyy = B_xX * D_xyy + B_xY * D_yyy + B_xZ * D_zyy - B_ * D_xyy;
      S.D_xyz = B_xX * D_xyz + B_xY * D_yyz + B_xZ * D_zyz - B_ * D_xyz;
      S.D_xzz = B_xX * D_xzz + B_xY * D_yzz + B_xZ * D_zzz - B_ * D_xzz;
      S.D_yxx = B_yX * D_xxx + B_yY * D_yxx + B_yZ * D_zxx - B_ * D_yxx;
      S.D_yxy = B_yX * D_xxy + B_yY * D_yxy + B_yZ * D_zxy - B_ * D_yxy;
      S.D_yxz = B_yX * D_xxz + B_yY * D_yxz + B_yZ * D_zxz - B_ * D_yxz;
      S.D_yyy = B_yX * D_xyy + B_yY * D_yyy + B_yZ * D_zyy - B_ * D_yyy;
      S.D_yyz = B_yX * D_xyz + B_yY * D_yyz + B_yZ * D_zyz - B_ * D_yyz;
      S.D_yzz = B_yX * D_xzz + B_yY * D_yzz + B_yZ * D_zzz - B_ * D_yzz;
      S.D_zxx = B_zX * D_xxx + B_zY * D_yxx + B_zZ * D_zxx - B_ * D_zxx;
      S.D_zxy = B_zX * D_xxy + B_zY * D_yxy + B_zZ * D_zxy - B_ * D_zxy;
      S.D_zxz = B_zX * D_xxz + B_zY * D_yxz + B_zZ * D_zxz - B_ * D_zxz;
      S.D_zyy = B_zX * D_xyy + B_zY * D_yyy + B_zZ * D_zyy - B_ * D_zyy;
      S.D_zyz = B_zX * D_xyz + B_zY * D_yyz + B_zZ * D_zyz - B_ * D_zyz;
      S.D_zzz = B_zX * D_xzz + B_zY * D_yzz + B_zZ * D_zzz - B_ * D_zzz;
      
      S.Theta_ = -B_ * Theta_ + T(0.5) * alpha_ * (T(2) * (A_x * (D_X - E_X - T(2) * Z_X) + A_y * (D_Y - E_Y - T(2) * Z_Y) + A_z * (D_Z - E_Z - T(2) * Z_Z)) + 
         D_xXX * Gamma_Xxx + D_xXY * Gamma_Xxy + D_xXZ * Gamma_Xxz + D_xXY * Gamma_Xxy + D_xYY * Gamma_Xyy + D_xYZ * Gamma_Xyz + D_xXZ * Gamma_Xxz + D_xYZ * Gamma_Xyz + D_xZZ * Gamma_Xzz + 
         D_yXX * Gamma_Yxx + D_yXY * Gamma_Yxy + D_yXZ * Gamma_Yxz + D_yXY * Gamma_Yxy + D_yYY * Gamma_Yyy + D_yYZ * Gamma_Yyz + D_yXZ * Gamma_Yxz + D_yYZ * Gamma_Yyz + D_yZZ * Gamma_Yzz + 
         D_zXX * Gamma_Zxx + D_zXY * Gamma_Zxy + D_zXZ * Gamma_Zxz + D_zXY * Gamma_Zxy + D_zYY * Gamma_Zyy + D_zYZ * Gamma_Zyz + D_zXZ * Gamma_Zxz + D_zYZ * Gamma_Zyz + D_zZZ * Gamma_Zzz - (
         K_Xx * K_Xx + K_Xy * K_Yx + K_Xz * K_Zx + K_Yx * K_Xy + K_Yy * K_Yy + K_Yz * K_Zy + K_Zx * K_Xz + K_Zy * K_Yz + K_Zz * K_Zz + 
         D_X * (D_x - T(2) * Z_x) + D_Y * (D_y - T(2) * Z_y) + D_Z * (D_z - T(2) * Z_z) + K_ * Q_)
      );
      
      S.Z_x = -B_ * Z_x + B_xX * Z_x + B_yX * Z_y + B_zX * Z_z + alpha_ * (A_x * Q_ + (D_x - A_x - T(2) * Z_x) * K_xx + (D_y - A_y - T(2) * Z_y) * K_xy + (D_z - A_z - T(2) * Z_z) * K_xz - 
         (K_xx * Gamma_Xxx + K_xy * Gamma_Yxx + K_xz * Gamma_Zxx + K_xy * Gamma_Xxy + K_yy * Gamma_Yxy + K_yz * Gamma_Zxy + K_xz * Gamma_Xxz + K_yz * Gamma_Yxz + K_zz * Gamma_Zxz)
      );
      S.Z_y = -B_ * Z_y + B_xY * Z_x + B_yY * Z_y + B_zY * Z_z + alpha_ * (A_y * Q_ + (D_x - A_x - T(2) * Z_x) * K_xy + (D_y - A_y - T(2) * Z_y) * K_yy + (D_z - A_z - T(2) * Z_z) * K_yz - 
         (K_xx * Gamma_Xxy + K_xy * Gamma_Yxy + K_xz * Gamma_Zxy + K_xy * Gamma_Xyy + K_yy * Gamma_Yyy + K_yz * Gamma_Zyy + K_xz * Gamma_Xyz + K_yz * Gamma_Yyz + K_zz * Gamma_Zyz)
      );
      S.Z_z = -B_ * Z_z + B_xZ * Z_x + B_yZ * Z_y + B_zZ * Z_z + alpha_ * (A_z * Q_ + (D_x - A_x - T(2) * Z_x) * K_xz + (D_y - A_y - T(2) * Z_y) * K_yz + (D_z - A_z - T(2) * Z_z) * K_zz - 
         (K_xx * Gamma_Xxz + K_xy * Gamma_Yxz + K_xz * Gamma_Zxz + K_xy * Gamma_Xyz + K_yy * Gamma_Yyz + K_yz * Gamma_Zyz + K_xz * Gamma_Xzz + K_yz * Gamma_Yzz + K_zz * Gamma_Zzz)
      );
      
      return S;
   }
   
private:

   std::array<T, 50> U;
   
   T& alpha_;
   T beta_X;
   T beta_Y;
   T beta_Z;
   
   T& gamma_xx;
   T& gamma_xy;
   T& gamma_xz;
   T& gamma_yy;
   T& gamma_yz;
   T& gamma_zz;
   
   T& K_xx;
   T& K_xy;
   T& K_xz;
   T& K_yy;
   T& K_yz;
   T& K_zz;
   
   T& A_x;
   T& A_y;
   T& A_z;
   
   T& B_xX;
   T& B_xY;
   T& B_xZ;
   T& B_yX;
   T& B_yY;
   T& B_yZ;
   T& B_zX;
   T& B_zY;
   T& B_zZ;
   
   T& D_xxx;
   T& D_xxy;
   T& D_xxz;
   T& D_xyy;
   T& D_xyz;
   T& D_xzz;
   T& D_yxx;
   T& D_yxy;
   T& D_yxz;
   T& D_yyy;
   T& D_yyz;
   T& D_yzz;
   T& D_zxx;
   T& D_zxy;
   T& D_zxz;
   T& D_zyy;
   T& D_zyz;
   T& D_zzz;
   
   T& Theta_;
   
   T& Z_x;
   T& Z_y;
   T& Z_z;
   
};
