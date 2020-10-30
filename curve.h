//Известные параметры
#define a_pr    "-2835"
#define b_pr    "9774"
#define p_pr    "115792089237316195423570985008687907853269984665640564039457584007913111864739"
//посчитано на Python
#define d_pr    "3"
#define x_pr    "8119021769418921565894388378008753471190077835025535568042996557746564092404"
#define y_pr    "28771053522948763457620123068270166170394398578116419877580784537256379313025"
#define u_pr    "93528601524582384978654134272795468222405403055717890280271688132874849008326"
#define v_pr    "14443324612566128911211262381388707474030458136470034119105598903952521080679"
#define q_pr    "115792089237316195423570985008687907852907080286716537199505774922796921406320"




//структура точки
struct  Point
{
     mpz_t x_cor;
     mpz_t y_cor;
     mpz_t z_cor;
};
//структура для параметров
struct  Params
{
     mpz_t p;
     mpz_t Q;
     mpz_t D;
};


//инициализация точки
void point_initialization(struct Point *dot,mpz_t x, mpz_t y, mpz_t z);
//суммирование точек
void point_sum(struct Point dot_1, struct Point dot_2,struct Point dot_3);
//удвоение точек
void point_double(struct Point dot_1,struct Point dot_3);
//проверка точка на прямой или нет (возвращает 0 если на кривой);
int CheckPoint(struct Point P);
//инициализация параметров
void params_initialization(struct Params *par);
//лесенка Монтгомери
struct Point ladder(mpz_t k,struct Point P);
