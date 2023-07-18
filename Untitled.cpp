#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <cmath>
using namespace std;
//*****************//
//вспомогательные вещи
vector<double> const operator - (const vector<double>& Rhs, const vector<double>& Lhs){
    vector<double> res(Rhs.size());
    for(int i=0; i< Rhs.size(); i++){
        res[i] = Rhs[i] - Lhs[i];
    }
    return res;
}

vector<double> const operator + (const vector<double>& Rhs, const vector<double>& Lhs){
    vector<double> res(Rhs.size());
    for(int i=0; i< Rhs.size(); i++){
        res[i] = Rhs[i] + Lhs[i];
    }
    return res;
}

vector<double> const operator / (const vector<double>& Lhs, double num){
    vector<double> res(Lhs.size());
    for(int i=0; i< Lhs.size(); i++){
        res[i] = Lhs[i] / num;
    }
    return res;
}

vector<double> const operator *  (const vector<double>& Lhs, double num){
    vector<double> res(Lhs.size());
    for(int i=0; i< Lhs.size(); i++){
        res[i] = Lhs[i] * num;
    }
    return res;
}
//dot product of vectors
double const operator ^ (const vector<double>& Rhs, const vector<double>& Lhs){
    double res=0;
    for(int i=0; i< Rhs.size(); i++){
        res += Rhs[i]*Lhs[i];
    }
    return res;
}

//*****************//

//генератор случ вещественных чисел
// в промежутке от -0.5 до 0.5 ==> A = sqrt(12)
double DoubleRand() {
    double upper = 0.5;
    double lower = -0.5;
    return (upper - lower) * ( (double)rand() / (double)RAND_MAX ) + lower;
};

//может не совсем хорошо использовать мощную структуру -- вектор,
//но попробуем так 
struct Atom {
    vector<double> positionOld ={0,0,0}; // X_Y_Z
    vector<double> position ={0,0,0}; // X_Y_Z
    vector<double> positionReal ={0,0,0}; // X_Y_Z вычисляется через скорость без учета периодичности
    vector<double> speed = {0,0,0}; // X_Y_Z
    vector<double> acceleration = {0,0,0}; // X_Y_Z
};

vector<double> calculateAverVelocityAxes(vector<Atom>& atoms, int atomNumber){
    vector<double> averVelocity = {0,0,0};
    for(auto& atom : atoms){
        averVelocity = averVelocity + atom.speed/atomNumber;
    }
    //cout<<averVelocity[0]<<" "<<averVelocity[1]<<" "<<averVelocity[2]<<endl;
    return averVelocity;
}

void setInitialCond (vector<Atom>& atoms, double max_X, double max_Y, double max_Z, double dt ) {
    // равномерно распределим шарки по объему c расстоян-м = 2
    // 10x10x10 => 125 частичек
    // захардкодим параметры
    int i=0;
    for(double x =1.0; x <= max_X-1; x+=2){
        for(double y =1.0; y <= max_Y-1; y+=2){
            for(double z =1.0; z<= max_Z-1; z+=2){
                atoms[i].position = {x, y, z};
                atoms[i].positionReal = {x, y, z};
                i++;
            }
        }
    }

    //рандомно раскидываем скорости {vx,vy,vz}*A; 
    
    for(auto& atom : atoms){
        atom.speed = {DoubleRand(), DoubleRand(), DoubleRand()};
    }
    

    //устраняем движение центра масс
    vector<double> averVelocity = calculateAverVelocityAxes( atoms ,125.);
    //A = sqrt(12) т.к. промежуток (-0.5 ; 0.5) 
    for(auto& atom : atoms){
        double factor = sqrt(12);
        atom.speed = (atom.speed - averVelocity)*factor; //перемасштабирование скорости
    }
    //алгоритм Верле требует прошлое значение позиции
    //зададим положения в 0-dt момент времени
    for(auto& atom : atoms) {
        atom.positionOld = atom.position - atom.speed * dt;
    }

}

vector<double> calculateRadiusVect(Atom& A, Atom& B) {
    //вычисление радиус вектора двух частиц с учетом периодич. гран условий
    //параметры ящика
    double max_X = 10;
    double max_Y = 10;
    double max_Z = 10;
    //затравочный рад-вектор
    vector<double> r0 = A.position - B.position;
    //учет периодичности
    r0[0] = r0[0] - max_X * round(r0[0] / max_X);
    r0[1] = r0[1] - max_Y * round(r0[1] / max_Y);
    r0[2] = r0[2] - max_Z * round(r0[2] / max_Z);
    return r0;
}

double radiusModule(vector<double>& r){
    return sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
}
double squareRadiusModule(vector<double>& r){
    return r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
}
double squareVelocityModule(vector<double>& r){
    return r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
}

void setAccelerationsToZero(vector<Atom>& atoms) {
    for(auto& atom : atoms){
        atom.acceleration = {0,0,0};
    }
}

//потенциальную энергию передадим по ссылке, чтоб изменять настоящее значение в main()
void calculateForcesAndPEnergy ( vector<Atom>& atoms, int N, double& currentPEnergy, uint64_t& collisionCount ) {
    //сбрасываем ускорения
    setAccelerationsToZero( atoms );
    //сбрасываем текущее значение потенциальной энергии
    currentPEnergy = 0;
    //параметр обрезки
    double rc = 3;
    //критическая сила по модулю
    double criticF_mod = (2/pow(rc,6) - 1) * 24/pow(rc, 7);
    //критический потенциал
    double criticV = (1/pow(rc,6) - 1) * 4/pow(rc,6);
    //перебор соседей
    for( int i=0; i< N-1; i++){
        for(int j=i+1; j< N; j++){
            vector<double> r = calculateRadiusVect(atoms[i], atoms[j]);
            double mod_r = radiusModule(r);
            //учитыываем соседей внутри сферы rc
            if(mod_r <= rc) {
                double pureF_mod = (24/pow(mod_r,7)) * (2/pow(mod_r,6) - 1); //Леннард Джонс сила
                double modifidedF_mod = pureF_mod - criticF_mod;
                // m=1 ==> F = a
                // f_ij = - f_ji
                atoms[i].acceleration = atoms[i].acceleration + (r/mod_r) * modifidedF_mod;
                atoms[j].acceleration = atoms[j].acceleration - (r/mod_r) * modifidedF_mod;
                //вычисляем потенциальную энергию
                currentPEnergy = currentPEnergy + (1/pow(mod_r,6) - 1) * 4/pow(mod_r,6) - criticV;
                //сигма = 1 ==> соударение при радиусе 1
                if(mod_r < 1) { collisionCount += 1; }
                
            }
        }
    }
}

void updatePositions (vector<Atom>& atoms, double dt) {
    //по алгоритму Верле r(t+dt) = 2r(t) - r(t-dt) +f/m *dt^2
    //параметры ящика
    double max_X = 10;
    double max_Y = 10;
    double max_Z = 10;
    //не забыть про периодические гран условия
    for(auto& atom : atoms) {
        vector<double> tmpCurrentPosition = atom.position;
        atom.position = atom.position * 2 - atom.positionOld + atom.acceleration *dt*dt;
        atom.positionOld = tmpCurrentPosition;
        //учет периодичности
        atom.position[0] -= max_X * floor(atom.position[0] / max_X);
        atom.position[1] -= max_Y * floor(atom.position[1] / max_Y);
        atom.position[2] -= max_Z * floor(atom.position[2] / max_Z);
        //обновление положений positionReal -- через скорости без перииодичности
        atom.positionReal = atom.positionReal + atom.speed*dt;
    }
    
}

void updateVelocity(vector<Atom>& atoms, vector<vector<double>>& old_acc ,double dt, int N) {
    //v(t+dt) = v(t) + (f(t-dt) + f(t))/2 *dt верле.
    for(int i=0; i<N; i++) {
        atoms[i].speed = atoms[i].speed + (old_acc[i] + atoms[i].acceleration)/2 *dt;
    }
}

double calculateMSD(vector<Atom>& atoms, vector<vector<double>>& initPositions) {
    double MSD = 0;
    for(int i =0; i< atoms.size(); i++) {
        vector<double> dr = atoms[i].positionReal - initPositions[i];
        MSD += squareRadiusModule(dr);
    }
    MSD = MSD/atoms.size();
    return MSD;
}

double calculateAVCF(vector<Atom>& atoms, vector<vector<double>>& initVelocities) {
    // ^ operator --  dot product of vectors
    double AV = 0;
    for(int i=0; i< atoms.size(); i++) {
        AV += atoms[i].speed ^ initVelocities[i];
    }
    AV = AV/atoms.size();
    return AV;
}

void calculateKinetEnergy(vector<Atom>& atoms, double& currentKEnergy) {
    //сбрасываем значение кинетиеской энергии
    currentKEnergy = 0;
    //вычисляем новое значение
    for(auto& atom : atoms) {
        double v2  = squareVelocityModule(atom.speed);
        currentKEnergy += v2 / 2;
    }

}

map<int, int> maxwellDistribution(vector<Atom>& atoms) {
    map<int, int> distribution; // уничтожится при выходе из функции
    for(auto& atom : atoms) {
        double v2 = squareVelocityModule(atom.speed);
        distribution[round(v2)] +=1;
    }
  return distribution;
}


//currentAVCF и currentMSD передаем по ссылке т.к. хотим менять глобальные значения
void makeTimeStep(vector<Atom>& atoms, double dt, int N, 
                    double& currentAVCF, double& currentMSD, vector<vector<double>>& initVelocities, vector<vector<double>>& initPositions,
                    double& currentPEnergy, double& currentKEnergy, uint64_t& collisionCount) {
    // v(t+dt) = v(t) + [a(t) + a(t+dt)]/2 *dt
    // сохраняем старые ускорения
    vector<vector<double>> old_acc(N);
    for(int i=0; i<N; i++) {
        old_acc[i] = atoms[i].acceleration;
    }
    //обновляем местоположения (те, которые без учета периодичности тоже)
    updatePositions (atoms, dt);
    currentMSD = calculateMSD( atoms, initPositions );//MSD вычисляется через positionReal
    //обновляем ускорения (силы) и вычисляем энергии
    calculateForcesAndPEnergy ( atoms, N, currentPEnergy, collisionCount);
    calculateKinetEnergy ( atoms, currentKEnergy );
    //обновляем скорости
    updateVelocity( atoms,  old_acc, dt,  N);
    currentAVCF = calculateAVCF( atoms, initVelocities );
    calculateAverVelocityAxes( atoms, N ); //проверил ЗСИ
}

int main(){
    //шаг времени
    double dt = 0.005;
    double t=0;
    //количество итераций
    int maxIterations = 10000;
    //параметры ящика
    double max_X = 10;
    double max_Y = 10;
    double max_Z = 10;
    //количество частиц
    int N = 125; 
    //количество соударений
    uint64_t collisionCount = 0;
    //массив всех частц
    vector<Atom> atoms(N);
    //мы не будем сохранять значения AVCF и MSD в массивы
    //организуем сразу записи в файл
    double currentAVCF = 0;
    double currentMSD = 0;
    //то же самое с энергиями
    double currentKEnergy = 0;
    double currentPEnergy = 0;
    //
    //поток для записи в файл
    ofstream out;         //поток для позиций
    ofstream outMSD;      //поток для MSD
    ofstream outAVCF;     //поток для AVCF
    ofstream outENERGY;     //поток для вывода энергий
    ofstream outMaxwell;    //поток для распределения скоростей
    outMaxwell.open("Maxwell.txt");
    outENERGY.open("Energies.txt");
    outAVCF.open("AVCF.txt");
    outMSD.open("MSD.txt");
    out.open("data.txt");
    if (out.is_open())
    {
        setInitialCond ( atoms, max_X, max_Y, max_Z, dt );
        //сохрняем начальные позиции и скорости для MSD  и AVCF
        //и выводим в файл начальные позиции
        //все одним циклом
        vector<vector<double>> initVelocities(N);
        vector<vector<double>> initPositions(N);
        for(int i=0; i<N; i++){
            initPositions[i] = atoms[i].position;
            initVelocities[i] = atoms[i].speed;
            out << ' '<< atoms[i].position[0]<<' '<< atoms[i].position[1]<<' '<< atoms[i].position[2]<< endl;
        }
        out<<endl;
        out<<endl;
        for(int j=0; j<maxIterations; j++){
            cout<<j<<endl;
            // основная функция
            makeTimeStep( atoms,  dt,  N, currentAVCF, currentMSD, initVelocities, initPositions, currentPEnergy, currentKEnergy, collisionCount);
            double fullEnergy = currentPEnergy + currentKEnergy; //полная энергия
            map<int,int> distribution = maxwellDistribution( atoms ); //распределение скоростей
            //записываем в файлы
            t+=dt;
            if(outMaxwell.is_open()) {
                for(auto& item : distribution){
                    outMaxwell<<item.first<<' '<<item.second<<endl;
                }
                outMaxwell<<endl;
                outMaxwell<<endl;
            }
            if (outENERGY.is_open()) {
                outENERGY<<fullEnergy<<" "<<currentPEnergy<<" "<<currentKEnergy<<" "<<t<<endl;
            }
            if (outMSD.is_open()) {
                outMSD<<currentMSD<<" "<<t<<endl;
            }
            if (outAVCF.is_open()) {
                outAVCF<<currentAVCF<<" "<<t<<endl;
            }

            for(int i=0; i<N; i++){
                out <<' '<<  atoms[i].position[0]<<' '<< atoms[i].position[1]<<' '<< atoms[i].position[2]<< endl;
            }
            out<<endl;
            out<<endl;
        }
    }
    cout<<"time: "<<(double)collisionCount/maxIterations;
    outMaxwell.close();
    outENERGY.close();
    outMSD.close();
    outAVCF.close();
    out.close();        

    return 0;
}