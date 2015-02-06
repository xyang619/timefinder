
#include <iostream>
#include <sstream>

using namespace std;

class A{
public:
    virtual ~A(){};
    virtual int f() const = 0;
};

class B: public A
{
public:
    B(int x);
    virtual ~B(){};
    int getX() const;
    void setX(int x);
private:
    int x;
};

class C: public B
{
public:
    C(int x);
    ~C(){};
    int f() const;
};

class C1: public B
{
public:
    C1(int x);
    ~C1(){};
    int f() const;
};

B::B(int x):x(x){ }

int B::getX() const{ return x;}

void B::setX( int x){
    this->x=x;
}

C::C(int x):B(x){ }

C1::C1(int x):B(x){ }

int C::f() const{
    return getX();
}

int C1::f() const{
    return getX()*getX();
}

int main(){
    B *t1 = new C(1);
    A *t2 = new C1(2);
    cout << t1->f() <<" " << t2->f();
    t1->setX(10);
    cout << t1->f() <<endl;
    delete t1;
    delete t2;
}