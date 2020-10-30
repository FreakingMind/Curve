#include <gmp.h>
#include "curve.h"


//инициализация точки
void point_initialization(struct Point *dot,mpz_t x, mpz_t y, mpz_t z){
    mpz_init_set_ui(dot->x_cor,x);
    mpz_init_set_ui(dot->y_cor,y);
    mpz_init_set_ui(dot->z_cor,z);

}
//инициализация параметров
void params_initialization(struct Params *par){
//    mpz_inits(par->U,par->V,pa    r->Q,par->D);
    mpz_init_set_str(par->p,p_pr,10);
    mpz_init_set_str(par->Q,q_pr,10);
    mpz_init_set_str(par->D,d_pr,10);
}
//суммирование точек A = X1Y2, B = Y1X2,X3 = BY1-Y2A, Y3 = X1A-BX2, Z3 = Y2X2-X1Y1
void point_sum(struct Point dot_1, struct Point dot_2,struct Point dot_3){
    mpz_t A;
    mpz_t B;
    mpz_t X3;
    mpz_t Y3;
    mpz_t Z3;
    mpz_t prob;
    mpz_init_set_str(prob,p_pr,10);
    mpz_init(A);
    mpz_init(B);
    mpz_init(X3);
    mpz_init(Y3);
    mpz_init(Z3);

    mpz_mul(A,dot_1.x_cor,dot_2.y_cor);
    mpz_mul(B,dot_1.y_cor,dot_2.x_cor);

    //x3
    mpz_mul(B,B,dot_1.y_cor);
    mpz_mul(A,A,dot_2.y_cor);
    mpz_sub(X3,B,A);
    mpz_init(A);
    mpz_init(B);
    //y3
    mpz_mul(A,dot_1.x_cor,dot_2.y_cor);
    mpz_mul(B,dot_1.y_cor,dot_2.x_cor);
    mpz_mul(A,A,dot_1.x_cor);
    mpz_mul(B,B,dot_2.x_cor);
    mpz_sub(Y3,A,B);

    mpz_init(A);

    mpz_init(B);
    //z3
    mpz_mul(A,dot_2.y_cor,dot_2.x_cor);
    mpz_mul(B,dot_1.x_cor,dot_1.y_cor);
    mpz_sub(Z3,A,B);
    mpz_mod(X3,X3,prob);
    mpz_mod(Y3,Y3,prob);
    mpz_mod(Z3,Z3,prob);


//    mpz_set(dot_3.x_cor,X3);
//    mpz_set(dot_3.y_cor,Y3);
//    mpz_set(dot_3.z_cor,Z3);
    point_initialization(&dot_3,X3,Y3,Z3);
    mpz_clear(A);
    mpz_clear(B);
    mpz_clear(X3);
    mpz_clear(Y3);
    mpz_clear(Z3);

}

//удвоение точки
//A=X1^2, B=Y1^2, D=A+B, G=(X1+Y1)^2-D, X3=(2Y1-G)x(X1+A+1), Y3 = (G-2X1)x(Y1+B+1), Z3 = (X1-Y1)x(G+2D)
void point_double(struct Point dot_1,struct Point dot_3){
    mpz_t A;
    mpz_t B;
    mpz_t D;
    mpz_t G;
    mpz_t X3;
    mpz_t Y3;
    mpz_t Z3;
    mpz_t bonus;
    mpz_t prob;
    mpz_init_set_str(prob,p_pr,10);
    mpz_init(A);
    mpz_init(B);
    mpz_init(D);
    mpz_init(G);
    mpz_init(X3);
    mpz_init(Y3);
    mpz_init(Z3);
    mpz_init(bonus);

    //A
    mpz_mul(A,dot_1.x_cor,dot_1.x_cor);
    //B
    mpz_mul(B,dot_1.y_cor,dot_1.y_cor);
    //G
    mpz_add(G,dot_1.x_cor,dot_1.y_cor);
    mpz_mul(G,G,G);
    mpz_sub(G,G,D);
    //X3
    //mpz_mul_2exp(X3,dot_1.y_cor,1);
    mpz_mul(X3,dot_1.y_cor,2);
    mpz_sub(X3,X3,G);
    mpz_add(Y3,A,1);
    mpz_add(Y3,Y3,dot_1.x_cor);
    mpz_mul(X3,X3,Y3);

    mpz_init(Y3);
    //Y3
    mpz_mul(bonus,dot_1.x_cor,2);
    mpz_sub(bonus,G,bonus);
    mpz_add(Y3,dot_1.y_cor,B);
    mpz_add(Y3,Y3,1);
    mpz_mul(Y3,bonus,Y3);

    mpz_init(bonus);
    //Z3
    mpz_sub(bonus,dot_1.x_cor,dot_1.y_cor);
    mpz_mul(Z3,D,2);
    mpz_add(Z3,Z3,G);
    mpz_mul(Z3,Z3,bonus);
    mpz_mod(X3,X3,prob);
    mpz_mod(Y3,Y3,prob);
    mpz_mod(Z3,Y3,prob);
    point_initialization(&dot_3,X3,Y3,Z3);
    mpz_clear(A);
    mpz_clear(B);
    mpz_clear(D);
    mpz_clear(G);
    mpz_clear(X3);
    mpz_clear(Y3);
    mpz_clear(Z3);
    mpz_clear(bonus);
}


//Лесенка Монтгомери
struct Point ladder(mpz_t k,struct Point P)
{   struct Point R,Q;
    point_initialization(&R,P.x_cor,P.y_cor,P.z_cor);
    point_initialization(&Q,1,-1,0);
    int bits_num;
    bits_num = mpz_sizeinbase(k,2);
    for (int i = bits_num-1;i>=0;i--)
    { if (mpz_tstbit(k,i) == 0){
            point_sum(R,Q,R);
            point_double(Q,Q);
        }
       else
        {
            point_sum(Q,R,Q);
            point_double(R,R);
        }
    }
    return Q;
}

//Проверка точки на прямой
//x^3+y^3+1=3*x*y
int CheckPoint(struct Point P){
    mpz_t left;
    mpz_t left_more;
    mpz_t right;
    mpz_t probik;
    mpz_init_set_str(probik,p_pr,10);
    mpz_init(left);
    mpz_init(right);
    mpz_init(left_more);

    mpz_add(left,left,P.x_cor);
    mpz_mul(right,left,left);
    mpz_mul(left,right,left); //x^3
    mpz_init(right);
    mpz_add(left_more,left_more,P.y_cor);
    mpz_mul(right,left_more,left_more);
    mpz_mul(left_more,right,left_more); //y^3
    mpz_add(left,left,left_more);
    mpz_add(left,left,1); //x^3+y^3+1
    mpz_mod(left,left,probik);
    mpz_init(right);

    mpz_mul(right,P.y_cor,P.x_cor);
    mpz_mul(right,right,3);
    mpz_mod(right,right,probik);
    int check = mpz_cmp(right,left);
    if (check == 0)
    {
        mpz_clear(left);
        mpz_clear(probik);
        mpz_clear(right);
        mpz_clear(left_more);

    }
    return check;
}





















