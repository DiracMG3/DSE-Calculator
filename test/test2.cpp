
#include "init.h"

using namespace std;
using namespace Eigen;


int main()
{
    //clock_t start = clock();
    intfunc = f;
    intfunc1 = fa;
    
    for (int i = 0; i < n; i++)
    {
        A0(i) = 10;
        XX(i) = 0;
    }

   

    
    
    for (int j = 0; j < 100; j++)
    {
        FXold = FFX(A0);
        cout << "norm of fx is :    " << FXold.norm() << endl;

        jac = JacFX(A0);
        cout << "det of jac is :    " << jac.determinant() << endl;
        
        VX = jac.inverse() * FXold; 
        A1 = A0 - VX;
        cout << "iteration of:" << j << endl;
        FXnew = FFX(A1);
        cout << "norm of new fx is :    " << FXnew.norm() << endl;
        if (VX.norm() < 1e-8)
            break;
        
        cout << "norm of new VX is :    " << VX.norm() << endl;

        
        if (FXnew.norm() >= FXold.norm())
        {
            for (int i = 1; i < 100; i++)
            {
                A1 = A0 - pow(0.5,i) * VX;
                FXnew = FFX(A1);
                cout << setprecision(20) << FXold.norm() << "   " << FXnew.norm() << endl;
                if (FXnew.norm() < FXold.norm())
                    break;
            }
            A0 = A1;
        }
        else 
        {
            A0 = A1;
        }
        
    }
    
    ofstream out1("jac.dat");
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            out1 << jac(i,j) << "\t" << endl;
        }
        
    }
    out1.close();


    ofstream out("XX.dat");
    for (int i = 0; i < n; i++)
    {
        out << XX(i) << "\t" << endl;
    }
    out.close();

    // clock_t finish = clock();
    // double TotalTime = (double)(finish - start) / CLOCKS_PER_SEC;
    // cout << "\n程序的运行时间为" << TotalTime << "秒！" << endl;

    return 0;
}









