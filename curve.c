#include "curve.h"
#include <string.h>
#include <stdio.h>

//инициализация точки
void point_initialization(struct Point *dot,mpz_t x, mpz_t y, mpz_t z){
    mpz_init_set(dot->x_cor,x);
    mpz_init_set(dot->y_cor,y);
    mpz_init_set(dot->z_cor,z);
}

//инициализация параметров
void params_initialization(struct Params *par){
    mpz_init_set_str(par->p,p_pr,10);
    mpz_init_set_str(par->Q,q_pr,10);
    mpz_init_set_str(par->D,d_pr,10);
}

//суммирование точек
//X3 = Y1X2*Y1Z2-Z1Y2*X1Y2
//Y3 = X1Z2*X1Y2-Y1X2*Z1X2
//Z3 = Z1Y2*Z1X2-X1Z2*Y1Z2
struct Point point_sum(struct Point dot_1, struct Point dot_2,mpz_t P) {
    struct Point dot_3;
    mpz_t help;
    mpz_init(dot_3.x_cor);
    mpz_init(dot_3.y_cor);
    mpz_init(dot_3.z_cor);
    mpz_init(help);

    //x3
    mpz_mul(help,dot_1.y_cor,dot_2.x_cor);
    mpz_mul(help,help,dot_1.y_cor);
    mpz_mul(dot_3.x_cor,help,dot_2.z_cor);
    mpz_mul(help,dot_1.z_cor,dot_2.y_cor);
    mpz_mul(help,help,dot_2.y_cor);
    mpz_mul(help,help,dot_1.x_cor);
    mpz_sub(dot_3.x_cor,dot_3.x_cor,help);

    //y3
    mpz_mul(help,dot_1.x_cor,dot_2.z_cor);
    mpz_mul(help,help,dot_1.x_cor);
    mpz_mul(dot_3.y_cor,help,dot_2.y_cor);
    mpz_mul(help,dot_1.y_cor,dot_2.x_cor);
    mpz_mul(help,help,dot_1.z_cor);
    mpz_mul(help,help,dot_2.x_cor);
    mpz_sub(dot_3.y_cor,dot_3.y_cor,help);

    //z3
    mpz_mul(help,dot_1.z_cor,dot_2.y_cor);
    mpz_mul(help,help,dot_1.z_cor);
    mpz_mul(dot_3.z_cor,help,dot_2.x_cor);
    mpz_mul(help,dot_1.x_cor,dot_2.z_cor);
    mpz_mul(help,help,dot_2.z_cor);
    mpz_mul(help,help,dot_1.y_cor);
    mpz_sub(dot_3.z_cor,dot_3.z_cor,help);

    mpz_fdiv_r(dot_3.x_cor, dot_3.x_cor, P);
    mpz_fdiv_r(dot_3.y_cor, dot_3.y_cor, P);
    mpz_fdiv_r(dot_3.z_cor, dot_3.z_cor, P);
    mpz_clear(help);
    return(dot_3);
}

//удвоение точки
//X3 = Y1*(Z1^3-X1^3)
//Y3 = X1*(Y1^3-Z1^3)
//Z3 = Z1*(X1^3-Y1^3)
struct Point point_double(struct Point dot_1,mpz_t P){
    struct Point dot_2;
    mpz_t help;
    mpz_init(dot_2.x_cor);
    mpz_init(dot_2.y_cor);
    mpz_init(dot_2.z_cor);
    mpz_init(help);

    //x3
    mpz_mul(help,dot_1.z_cor,dot_1.z_cor);
    mpz_mul(dot_2.x_cor,help,dot_1.z_cor);
    mpz_mul(help,dot_1.x_cor,dot_1.x_cor);
    mpz_mul(help,help,dot_1.x_cor);
    mpz_sub(dot_2.x_cor,dot_2.x_cor,help);
    mpz_mul(dot_2.x_cor,dot_2.x_cor,dot_1.y_cor);

    //y3
    mpz_mul(help,dot_1.y_cor,dot_1.y_cor);
    mpz_mul(dot_2.y_cor,help,dot_1.y_cor);
    mpz_mul(help,dot_1.z_cor,dot_1.z_cor);
    mpz_mul(help,help,dot_1.z_cor);
    mpz_sub(dot_2.y_cor,dot_2.y_cor,help);
    mpz_mul(dot_2.y_cor,dot_2.y_cor,dot_1.x_cor);

    //z3
    mpz_mul(help,dot_1.x_cor,dot_1.x_cor);
    mpz_mul(dot_2.z_cor,help,dot_1.x_cor);
    mpz_mul(help,dot_1.y_cor,dot_1.y_cor);
    mpz_mul(help,help,dot_1.y_cor);
    mpz_sub(dot_2.z_cor,dot_2.z_cor,help);
    mpz_mul(dot_2.z_cor,dot_2.z_cor,dot_1.z_cor);

    mpz_fdiv_r(dot_2.x_cor, dot_2.x_cor, P);
    mpz_fdiv_r(dot_2.y_cor, dot_2.y_cor, P);
    mpz_fdiv_r(dot_2.z_cor, dot_2.z_cor, P);
    mpz_clear(help);
    return(dot_2);
}


//Лесенка Монтгомери
struct Point ladder(mpz_t k,struct Point P) {
    struct Point R,Q;
    struct Params param;
    mpz_t zerox, zeroy, zeroz;
    mpz_init_set_si(zerox, 1);
    mpz_init_set_si(zeroy, -1);
    mpz_init_set_si(zeroz, 0);
    params_initialization(&param);
    point_initialization(&R,P.x_cor,P.y_cor,P.z_cor);
    point_initialization(&Q,zerox,zeroy,zeroz);
    int bits_num;
    bits_num = mpz_sizeinbase(k, 2);
    for (int i = bits_num-1;i>=0;i--)
    {
        if (mpz_tstbit(k,i) == 0){
            R = point_sum(R,Q, param.p);
            Q = point_double(Q,param.p);
        }
        else
        {
            Q = point_sum(Q,R,param.p);
            R = point_double(R,param.p);
        }
    }
    mpz_clear(zerox);
    mpz_clear(zeroy);
    mpz_clear(zeroz);
    return Q;
}

//Проверка точки на прямой
//x^3+y^3+z^3=3*D*x*y*z
int CheckPoint(struct Point P) {
    mpz_t left,left_more,right,probik,d,three;
    mpz_init_set_str(d,d_pr,10);
    mpz_init_set_ui(three, 3);
    mpz_init_set_str(probik,p_pr,10);
    mpz_init(left);
    mpz_init(right);
    mpz_init(left_more);

    mpz_mul(left,P.x_cor,P.x_cor);
    mpz_mul(left,left,P.x_cor);
    mpz_mul(left_more,P.y_cor,P.y_cor);
    mpz_mul(left_more,left_more,P.y_cor);
    mpz_add(left,left,left_more);
    mpz_mul(left_more,P.z_cor,P.z_cor);
    mpz_mul(left_more,left_more,P.z_cor);
    mpz_add(left,left,left_more);
    mpz_fdiv_r(left,left,probik);

    mpz_mul(right,three,P.x_cor);
    mpz_mul(right,right,P.y_cor);
    mpz_mul(right,right,P.z_cor);
    mpz_mul(right,right,d);
    mpz_fdiv_r(right,right,probik);

    gmp_printf("Левая часть, x^3+y^3+z^3:\n%Zd\n",left);
    gmp_printf("Правая часть, 3Dxyz:\n%Zd\n",right);
    int check = mpz_cmp(right,left);
    mpz_clear(left);
    mpz_clear(probik);
    mpz_clear(right);
    mpz_clear(left_more);
    return check;
}
