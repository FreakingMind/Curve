#include <stdio.h>
#include <gmp.h>
#include "curve.h"


void Test(){
//    printf("Параметры:");
//    printf("a:%s\n",a_pr);
//    printf("b:%s\n",b_pr);
//    printf("p:%s\n",p_pr);
//    printf("x:%s\n",x_pr);
//    printf("y:%s\n",y_pr);
//    printf("m:%s\n",q_pr);
//    printf("Параметры Гессе\n");
//    printf("D:%s\n",d_pr);
//    printf("x:%s\n",u_pr);
//    printf("y:%s\n",v_pr);
//    printf("z:1\n");
    //инициализация параметров
    mpz_t u1,v1,z1;
    mpz_init_set_str(u1,u_pr,10);
    mpz_init_set_str(v1,v_pr,10);
    mpz_init_set_ui(z1,1);

    //проверка принадлежности точки
    printf("Тест 1. Принадлежит ли точка кривой?\n");
    mpz_t k,n;
    mpz_init(k);
    mpz_init_set_str(n,"9999999999999999999999999999999999999999999999999999999999999999999999999999999999",10);
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    gmp_randseed_ui(state, 8473522);
    mpz_urandomm(k, state, n);
    struct Params Chekp;
    params_initialization(&Chekp);
    struct Point My_point;
    point_initialization(&My_point,u1,v1,z1);
    printf("Сгенерированное k:\n");
    gmp_printf("%Zd\n",k);
    if (CheckPoint(ladder(k,My_point)) == 0){
        printf("Да, точка на прямой\n");
    }
    else {
        printf("Нет, точка не на прямой\n");
    }
    mpz_clear(u1);
    mpz_clear(v1);
    mpz_clear(z1);
    mpz_clear(k);

    //Второй тест
    printf("Тест 2. Координаты полученной точки q[P]?\n");
    mpz_t qx;
    mpz_init_set_si(qx, 1);
    mpz_t increased_q;
    mpz_init_set_str(increased_q,q_pr,10);
    struct Point P;
    mpz_init(P.x_cor);
    mpz_init(P.y_cor);
    mpz_init(P.z_cor);
    P = ladder(increased_q,My_point);
    gmp_printf("x=%Zd\n",P.x_cor);
    gmp_printf("y=%Zd\n",P.y_cor);
    gmp_printf("z=%Zd\n",P.z_cor);
    mpz_clear(P.x_cor);
    mpz_clear(P.y_cor);
    mpz_clear(P.z_cor);

    //Третий тест
    printf("Тест 3. Выполняется ли [q+1]P=P и [q-1]P=-P?\n");
    mpz_add(increased_q,increased_q,qx);
    mpz_t decreased_q,a,b;
    mpz_inits(a,b);
    mpz_init_set_str(decreased_q,q_pr,10);
    mpz_sub(decreased_q,decreased_q,qx);

    struct Point reversed;
    point_initialization(&reversed,My_point.y_cor,My_point.x_cor,My_point.z_cor); //обратный элемент

    struct Point LeftCheck;
    mpz_init(LeftCheck.x_cor);
    mpz_init(LeftCheck.y_cor);
    mpz_init(LeftCheck.z_cor);
    LeftCheck = ladder(increased_q,My_point); //[q+1]P

    printf("%s\n", "Полученная точка равна точке P, если их координаты пропорциональны, т.к. работаем в проективном пространстве");
    printf("%s\n", "Первое равенство:");
    mpz_mul(a,My_point.x_cor,LeftCheck.y_cor);
    mpz_mul(b,My_point.y_cor,LeftCheck.x_cor);
    mpz_fdiv_r(a,a,Chekp.p);
    mpz_fdiv_r(b,b,Chekp.p);
    printf("x1*y2:\n");
    gmp_printf("%Zd\n",a);
    printf("y1*x2:\n");
    gmp_printf("%Zd\n",b);

    mpz_mul(a,My_point.x_cor,LeftCheck.z_cor);
    mpz_mul(b,My_point.z_cor,LeftCheck.x_cor);
    mpz_fdiv_r(a,a,Chekp.p);
    mpz_fdiv_r(b,b,Chekp.p);
    printf("x1*z2:\n");
    gmp_printf("%Zd\n",a);
    printf("z1*x2:\n");
    gmp_printf("%Zd\n",b);

    mpz_mul(a,My_point.y_cor,LeftCheck.z_cor);
    mpz_mul(b,My_point.z_cor,LeftCheck.y_cor);
    mpz_fdiv_r(a,a,Chekp.p);
    mpz_fdiv_r(b,b,Chekp.p);
    printf("y1*z2:\n");
    gmp_printf("%Zd\n",a);
    printf("z1*y2:\n");
    gmp_printf("%Zd\n",b);

    printf("%s\n", "Полученная точка равна точке -P=(y,x,z), если их координаты пропорциональны, т.к. работаем в проективном пространстве");
    printf("%s\n", "Второе равенство:");

    struct Point RightCheck;
    RightCheck = ladder(decreased_q,My_point); // [q-1]P

    mpz_mul(a,reversed.x_cor,RightCheck.y_cor);
    mpz_mul(b,reversed.y_cor,RightCheck.x_cor);
    mpz_fdiv_r(a,a,Chekp.p);
    mpz_fdiv_r(b,b,Chekp.p);
    printf("x1*y2:\n");
    gmp_printf("%Zd\n",a);
    printf("y1*x2:\n");
    gmp_printf("%Zd\n",b);

    mpz_mul(a,reversed.x_cor,RightCheck.z_cor);
    mpz_mul(b,reversed.z_cor,RightCheck.x_cor);
    mpz_fdiv_r(a,a,Chekp.p);
    mpz_fdiv_r(b,b,Chekp.p);
    printf("x1*z2:\n");
    gmp_printf("%Zd\n",a);
    printf("z1*x2:\n");
    gmp_printf("%Zd\n",b);

    mpz_mul(a,reversed.y_cor,RightCheck.z_cor);
    mpz_mul(b,reversed.z_cor,RightCheck.y_cor);
    mpz_fdiv_r(a,a,Chekp.p);
    mpz_fdiv_r(b,b,Chekp.p);
    printf("y1*z2:\n");
    gmp_printf("%Zd\n",a);
    printf("z1*y2:\n");
    gmp_printf("%Zd\n",b);

    mpz_clear(increased_q);
    mpz_clear(decreased_q);
    mpz_clear(qx);
    mpz_clear(reversed.x_cor);
    mpz_clear(reversed.y_cor);
    mpz_clear(reversed.z_cor);
    mpz_clear(RightCheck.x_cor);
    mpz_clear(RightCheck.y_cor);
    mpz_clear(RightCheck.z_cor);
    mpz_clear(LeftCheck.x_cor);
    mpz_clear(LeftCheck.y_cor);
    mpz_clear(LeftCheck.z_cor);

    //Четвертый тест
    printf("Тест 4. Выполняется ли: [k1+k2]P = [k1]P + [k2]P?\n");
    //Генерируем случайные значения к1 и к2
    mpz_t k1,k2,sum_k;
    mpz_init(k1);
    mpz_init(k2);
    mpz_init(sum_k);
    mpz_urandomm(k1, state, n);
    mpz_urandomm(k2, state, n);
    gmp_printf ("k1: %Zd\n", k1);
    gmp_printf ("k2: %Zd\n", k2);

    struct Point poin_k1;
    mpz_init(poin_k1.x_cor);
    mpz_init(poin_k1.y_cor);
    mpz_init(poin_k1.z_cor);
    poin_k1 = ladder(k1,My_point);

    struct Point poin_k2;
    mpz_init(poin_k2.x_cor);
    mpz_init(poin_k2.y_cor);
    mpz_init(poin_k2.z_cor);
    poin_k2 = ladder(k2,My_point);

    mpz_add(sum_k,k1,k2);

    struct Point sum_of_k;
    mpz_init(sum_of_k.x_cor);
    mpz_init(sum_of_k.y_cor);
    mpz_init(sum_of_k.z_cor);
    sum_of_k = ladder(sum_k,My_point);

    struct Point sums;
    mpz_init(sums.x_cor);
    mpz_init(sums.y_cor);
    mpz_init(sums.z_cor);
    sums = point_sum(poin_k1,poin_k2,Chekp.p);

    printf("Аналогично проверяем пропорциональность\n");
    printf("P[k1+k2].x * P[k1]+P[k2].y:\n");
    mpz_mul(a,sum_of_k.x_cor,sums.y_cor);
    mpz_fdiv_r(a,a,Chekp.p);
    gmp_printf("%Zd\n",a);
    printf("P[k1+k2].y * P[k1]+P[k2].x:\n");
    mpz_mul(b,sum_of_k.y_cor,sums.x_cor);
    mpz_fdiv_r(b,b,Chekp.p);
    gmp_printf("%Zd\n",b);

    printf("P[k1+k2].x * P[k1]+P[k2].z:\n");
    mpz_mul(a,sum_of_k.x_cor,sums.z_cor);
    mpz_fdiv_r(a,a,Chekp.p);
    gmp_printf("%Zd\n",a);
    printf("P[k1+k2].z * P[k1]+P[k2].x:\n");
    mpz_mul(b,sum_of_k.z_cor,sums.x_cor);
    mpz_fdiv_r(b,b,Chekp.p);
    gmp_printf("%Zd\n",b);

    printf("P[k1+k2].y * P[k1]+P[k2].z:\n");
    mpz_mul(a,sum_of_k.y_cor,sums.z_cor);
    mpz_fdiv_r(a,a,Chekp.p);
    gmp_printf("%Zd\n",a);
    printf("P[k1+k2].z * P[k1]+P[k2].y:\n");
    mpz_mul(b,sum_of_k.z_cor,sums.y_cor);
    mpz_fdiv_r(b,b,Chekp.p);
    gmp_printf("%Zd\n",b);

    mpz_clear(Chekp.p);
    mpz_clear(Chekp.D);
    mpz_clear(Chekp.Q);
    mpz_clear(My_point.x_cor);
    mpz_clear(My_point.y_cor);
    mpz_clear(My_point.z_cor);
    mpz_clear(k1);
    mpz_clear(k2);
    mpz_clear(sum_k);
    mpz_clear(n);
    mpz_clear(poin_k1.x_cor);
    mpz_clear(poin_k1.y_cor);
    mpz_clear(poin_k1.z_cor);
    mpz_clear(poin_k2.x_cor);
    mpz_clear(poin_k2.y_cor);
    mpz_clear(poin_k2.z_cor);
    mpz_clear(sum_of_k.x_cor);
    mpz_clear(sum_of_k.y_cor);
    mpz_clear(sum_of_k.z_cor);
    mpz_clear(sums.x_cor);
    mpz_clear(sums.y_cor);
    mpz_clear(sums.z_cor);
    gmp_randclear(state);
}


int main()
{
    Test();
    return 0;
}
