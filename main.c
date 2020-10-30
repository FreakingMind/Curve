#include <stdio.h>
#include <gmp.h>
#include "curve.h"


void Test(){
    printf("Параметры:");
    printf("a:%s\n",a_pr);
    printf("b:%s\n",b_pr);
    printf("p:%s\n",p_pr);
    printf("x:%s\n",x_pr);
    printf("y:%s\n",y_pr);
    printf("m:%s\n",q_pr);
    printf("Параметры Гессе\n");
    printf("D:%s\n",d_pr);
    printf("x:%s\n",u_pr);
    printf("y:%s\n",v_pr);
    printf("z:1\n");
//инициализация параметров
    mpz_t u1,v1,z1;
    mpz_init_set_str(u1,u_pr,10);
    mpz_init_set_str(v1,v_pr,10);
    mpz_init_set(z1,1);
// проверка принадлежности точки
//Генирируем случайную точку
    printf("Тест 1. Принадлежит ли точка кривой?\n");
    mpz_t k,n;
    mpz_init_set_ui(n,999999999999999999);
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    mpz_urandomm(k, state, n);
    struct Params Chekp;
    params_initialization(&Chekp);
    struct Point My_point;
    point_initialization(&My_point,u1,v1,z1);
    if (CheckPoint(ladder(k,My_point)) == 0){
        printf("Да, точка на прямой\n");
    }
    else{
        printf("Нет, точка не на прямой\n");
    }
//Второй тест
    printf("Тест 2. Равна ли нейтральная точка и q[P]?\n");
    struct Point qu;
    point_initialization(&qu,1,-1,0);
    mpz_t increased_q;
    mpz_init_set_str(increased_q,q_pr,10);
    struct Point P = ladder(increased_q,My_point);
    if ((P.x_cor == qu.x_cor) && (P.y_cor == qu.y_cor)&&(P.z_cor == qu.z_cor)){
        printf("Тест 2 пройден");
    }
    else{
        printf("Тест 2 не пройден");
    }

//Третий тест
    printf("Выполняется ли [q+1]P=P и [q-1]P = -P?");
    mpz_add(increased_q,increased_q,1);
    mpz_t decreased_q;
    mpz_init_set_str(decreased_q,q_pr,10);
    mpz_sub(decreased_q,decreased_q,1);

    struct Point reversed;
    point_initialization(&reversed,My_point.y_cor,My_point.x_cor,My_point.z_cor); //обратный элемент
    struct Point LeftCheck = ladder(increased_q,My_point); //[q+1]P
    struct Point RightCheck = ladder(decreased_q,My_point); // [q-1]P
    if ((LeftCheck.x_cor == My_point.x_cor)&& (RightCheck.y_cor == My_point.y_cor)&& (LeftCheck.z_cor == My_point.z_cor)&&
            (RightCheck.x_cor == reversed.x_cor) && (RightCheck.y_cor == reversed.y_cor) && (RightCheck.z_cor == reversed.z_cor))
    {
        printf("Тест 3 пройден.[q+1]P=P и [q-1]P = -P ");

    }
    else{
        printf("Тест 3 не пройден");
    }
//Четвертый тест
//Генерируем случайные точки к1 и к2
    mpz_t k1,k2,sum_k;
    mpz_init(sum_k);
    mpz_init_set_ui(n,99999999999999);
    gmp_randstate_t state_1;
    gmp_randinit_mt(state_1);
    gmp_randstate_t state_2;
    gmp_randinit_mt(state_2);
    mpz_urandomm(k1, state_1, n);
    mpz_urandomm(k2, state_2, n);

    struct Point poin_k1 = ladder(k1,My_point);
    struct Point poin_k2 = ladder(k2,My_point);
    mpz_add(sum_k,k1,k2);
    struct Point sum_of_k = ladder(sum_k,My_point);
    struct Point sums;
    point_sum(poin_k1,poin_k2,sums);
    if ((sum_of_k.x_cor == sums.x_cor)&&(sum_of_k.y_cor == sums.y_cor) && (sum_of_k.z_cor == sums.z_cor)){
        printf("Тест 4 пройден, [k1]P+k2[P] =[k1+k2]P");
    }
}


int main()
{
     Test();
    return 0;
}
