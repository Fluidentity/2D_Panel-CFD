#include <iostream>
using namespace std;



class Class_A {
private: 
    int a;
    int b;

public: 


Class_A() {
    cout << "Constructor of A" << endl;
}
Class_A(const Class_A &obj) {
    this->a = obj.a;
    this->b = obj.b;
}
Class_A(int x, int y) {this->a = x;this->b = y;}
virtual void func() {
        cout << "from Func_A" <<endl;
    }

};
class Class_B:public Class_A {
private:
    int c;
    int d;
public:
    Class_B() {
        cout << "Constrctor of B" << endl;
    }
    Class_B(int x, int y) {
        cout << "Constrctor of B" << endl;
    }
    void func() {
        cout << "from Func_B" << endl;
    }

};

int main() {
    int nx;
    float xMax = 0.5;
    float xMin = 0;
    float dx_temp = 0.2;


    nx = (int)((xMax-xMin)/dx_temp);

    cout << nx;





    return 0;

}