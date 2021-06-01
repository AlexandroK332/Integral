#include <iostream>
#include <string>
#include <random>
#include <ctime>

struct Res

{
    double Integ = 0;
    double pogr = 0;
   // Res operator=(Res a)
    //{
    //    Integ = a.Integ;
    //    pogr = a.pogr;

   // }
};
class Function

{
private:
    std::string F;
    std::string bufstr;
    double x;
    bool sravnznak(char a)
    
    {
        if ((a == '+') or (a == '*') or (a == '/') or (a == '^'))
        
        {
            return(0);
        
        }
        else
        
        {
            return(1);
        
        }
    
    }
    double Operation(double a, double b, char Z)
   
    {
        if (Z == '+')
        
        {
            return(a + b);
        
        }
        if (Z == '^')
        
        {
            return(pow(a, b));
        
        }
        if (Z == '*')
        
        {
            return(a * b);
        }
        if (Z == '/')
        
        {
            return(a / b);
        
        }
    }
    void Insert()
    {
        int a = F.length();
        for (size_t i = 0; i < a; i++)
        {
            if (F[i] == 'x')
            {
                F.erase(i, 1);
                F.insert(i, std::to_string(x));
                a = F.length();
            }
        }
    }
    void mFunction()
    {
        int a = F.length();
        for (size_t k = 0; k < a; k++)
        {
            int j = 0;
            for (size_t i = 0; i < a; i++)
            {
                if (F[i] == '(')
                {
                    j = i;
                   
                }
                if (F[i] == ')')
                {
                    Function q(F.substr(j + 1, i - j - 1), x);
                    F.erase(j, i - j + 1);
                    F.insert(j, std::to_string(q.GetRes(x)));
                    a = F.length();
                    break;
                }
              
            }
        }
    }
    void plus()
    {
        int a = F.length();
        for (size_t i = 0; i < a - 1; i++)
        {
            if ((F[i] == '+') && (F[i + 1] == '+'))
            {
                F.erase(i, 2);
                F.insert(i, "+");
                a = F.length();
            }
        }
    }
    void probel()
    {
        int a = F.length();
        for (size_t i = 0; i < a - 1; i++)
        {
            if ((F[i] == ' '))
            {
                F.erase(i, 1);
                a = F.length();
            }
        }

    }
    void binminus()
    {
        int a = F.length();
        for (size_t i = 0; i < a; i++)
        {
            if ((F[i] == '-') && (i > 0))
            {
                if ((F[i] == '-') && (sravnznak(F[i - 1])) && (F[i - 1] != '-')&&(F[i-1]!='s')&&(F[i-1]!='n'))
                {
                    F.insert(i, "+");
                    a = F.length();
                }
            }
        }
    }
    void twominus()
    {
        int a = F.length();
        for (size_t i = 0; i < a - 1; i++)
        {
            if ((F[i] == '-') && (F[i + 1] == '-'))
            {
                F.erase(i, 2);
                F.insert(i, "+");
                a = F.size();
            }
        }
    }
    void znak(char Z)
    {
        bool op0 = 1;
        int a = F.length();
        while (op0)
        {
            int count = 0;
            for (int i = 0; i < a; i++)
            {
                if (F[i] == Z)
                {
                    count++;
                    int buf1 = i - 1;
                    if (buf1 < 0)
                    {
                        throw 1;
                    }
                    while ((buf1 > 0) && (sravnznak(F[buf1])))
                    {
                        buf1--;
                    }
                    if (sravnznak(F[buf1]) == 0)
                    {
                        buf1++;
                    }
                    double o = std::stod(F.substr(buf1, i - buf1));
                    int buf2 = i + 1;
                    while (sravnznak(F[buf2]) && (buf2 != a))
                    {
                        buf2++;
                    }
                    double b = std::stod(F.substr(i + 1, buf2 - i + 1));
                    double c;
                    c = Operation(o, b, Z);
                    F.erase(buf1, buf2 - buf1);
                    F.insert(buf1, std::to_string(c));
                    a = F.size();
                }

            }
            if (count == 0)
            {
                op0 = 0;
            }
            a = F.length();
        }

    }
    void Trig()
    {
        int a = F.length();
        for (int i = 2; i < a; i++)
        {
            if ((F[i] == 'n') && (F[i - 1] == 'i') && (F[i - 2] == 's')&&(i+1<a))
            {
                int buf = i + 1;
                while (sravnznak(F[buf]) && (buf < a))
                {
                    buf++;
                }
                double b = sin(std::stod(F.substr(i+1, buf - i + 1)));
                F.erase(i - 2, buf - i + 1);
                F.insert(i - 2, std::to_string(b));
                a = F.length();
            }
            if ((F[i] == 's') && (F[i - 1] == 'o') && (F[i - 2] == 'c') && (i + 1 < a))
            {
                int buf = i + 1;
                while (sravnznak(F[buf]) && (buf < a))
                {
                    buf++;
                }
                double b = cos(std::stod(F.substr(i + 1, buf - i + 1)));
                F.erase(i - 2, buf - i+1);
                F.insert(i - 2, std::to_string(b));
                a = F.length();
            }
            if ((F[i] == 'n') && (F[i - 1] == 'a') && (F[i - 2] == 't') && (i + 1 < a))
            {
                int buf = i + 1;
                while (sravnznak(F[buf]) && (buf < a))
                {
                    buf++;
                }
                double b = tan(std::stod(F.substr(i + 1, buf - i + 1)));
                F.erase(i - 2, buf - i + 1);
                F.insert(i - 2, std::to_string(b));
                a = F.length();
            }
        }

    }
    
public:
    Function()
    {
        F = "x";
        x = 1;
        bufstr = F;
    }
    Function(std::string F_, double x_)
    {
        x = x_;
        F = F_;
        bufstr = F;
    }
    void SetStr(std::string str_)
    {
        F = str_;
        bufstr = F;
    }
    double GetRes(double x_)
    {

        x = x_;
        try
        {
            probel();
            Insert();
            mFunction();
            plus();
            binminus();
            twominus();
            Trig();
            znak('^');
            znak('*');
            znak('/');
            znak('+');
            double result = std::stod(F);
            F = bufstr;
            return(result);
        }
        catch (std::invalid_argument)
        {
            std::cerr << '\n' << "Input error" << std::endl;
            exit(1);
        }
        catch (int)
        {
            std::cerr << '\n' << "Input error" << std::endl;
            exit(1);
        }
    }
};
template <typename T>
class Integral
{
private:
    double toch;
    double h;
    T p1, p2;
    Res Res1;
    std::string str;
    Function A;
    bool sravnznak(char a)
    {
        if ((a == '+') or (a == '*') or (a == '/') or (a == '^'))
        {
            return(0);
        }
        else
        {
            return(1);
        }
    }
    double Formula(std::string F, double x)
    {
        int a = F.length();
        for (int i = 0; i < a; i++)
        {
            if (F[i] == 'x')
            {
                F.erase(i, 1);
                F.insert(i, std::to_string(x));
                a = F.length();
            }
        }
        a = F.length();
        for (int i = 0; i < a; i++)
        {
            if (F[i] == '(')
            {
                int j = i + 1;
                while ((F[j] != ')') && (j < a))
                {
                    j++;
                }
                double q = Formula(F.substr(i + 1, j - i - 1), x);
                F.erase(i, j - i + 1);
                F.insert(i, std::to_string(q));
                a = F.length();
            }
        }
        a = F.length();
        for (int i = 0; i < a; i++)
        {
            if ((F[i] == '-') && (i > 0))
            {
                if ((F[i] == '-') && (sravnznak(F[i - 1]) && (F[i - 1] != '-')))
                {
                    F.insert(i, "+");
                    a = F.size();
                }
            }
        }
        a = F.length();
        for (int i = 0; i < a - 1; i++)
        {
            if ((F[i] == '-') && (F[i + 1] == '-'))
            {
                F.erase(i, 2);
                F.insert(i, "+");
                a = F.size();
            }
        }
        for (int i = 0; i < a - 1; i++)
        {
            if ((F[i] == '+') && (F[i + 1] == '+'))
            {
                F.erase(i, 2);
                F.insert(i, "+");
                a = F.size();
            }
        }
        bool op0 = 1;
        while (op0)
        {
            int count = 0;
            for (int i = 0; i < a; i++)
            {
                if (F[i] == '^')
                {
                    count++;
                    int buf1 = i - 1;
                    while ((sravnznak(F[buf1])) && (buf1 != 0))
                    {
                        buf1--;
                    }
                    if (sravnznak(F[buf1]) == 0)
                    {
                        buf1++;
                    }
                    double o = std::stod(F.substr(buf1, i - buf1));
                    int buf2 = i + 1;
                    while (sravnznak(F[buf2]) && (buf2 != a))
                    {
                        buf2++;
                    }
                    double b = std::stod(F.substr(i + 1, buf2 - i + 1));
                    double c;
                    c = pow(o, b);


                    F.erase(buf1, buf2 - buf1);
                    F.insert(buf1, std::to_string(c));
                    a = F.size();
                }

            }
            if (count == 0)
            {
                op0 = 0;
            }
            a = F.length();
        }
        bool op1 = 1;
        while (op1)
        {
            int count = 0;
            for (int i = 0; i < a; i++)
            {
                if ((F[i] == '*') || (F[i] == '/'))
                {
                    count++;
                    int buf1 = i - 1;
                    while ((sravnznak(F[buf1])) && (buf1 != 0))
                    {
                        buf1--;
                    }
                    if (sravnznak(F[buf1]) == 0)
                    {
                        buf1++;
                    }
                    double o = std::stod(F.substr(buf1, i - buf1));
                    int buf2 = i + 1;
                    while (sravnznak(F[buf2]) && (buf2 != a))
                    {
                        buf2++;
                    }
                    double b = std::stod(F.substr(i + 1, buf2 - i + 1));
                    double c;
                    if (F[i] == '*')
                    {
                        c = o * b;
                    }
                    else
                    {

                        c = o / b;
                    }
                    F.erase(buf1, buf2 - buf1);
                    F.insert(buf1, std::to_string(c));
                    a = F.size();
                }

            }
            if (count == 0)
            {
                op1 = 0;
            }
            a = F.size();
        }
        bool op2 = 1;
        while (op2)
        {
            int count = 0;
            for (int i = 0; i < a; i++)
            {
                if ((F[i] == '+'))
                {
                    count++;
                    int buf1 = i - 1;
                    while ((sravnznak(F[buf1])) && (buf1 != 0))
                    {
                        buf1--;
                    }
                    if (sravnznak(F[buf1]) == 0)
                    {
                        buf1++;
                    }
                    double o = std::stod(F.substr(buf1, i - buf1));
                    int buf2 = i + 1;
                    while (sravnznak(F[buf2]) && (buf2 != a))
                    {
                        buf2++;
                    }
                    double b = std::stod(F.substr(i + 1, buf2 - i));
                    double c;

                    c = o + b;

                    F.erase(buf1, buf2 - buf1);
                    F.insert(buf1, std::to_string(c));
                    a = F.size();
                }

            }
            if (count == 0)
            {
                op2 = 0;
            }
            a = F.size();
        }
        return(std::stod(F));
    }
    double Gaussm(double p1, double p2)
    {
        const double Xi[3] = { -0.7745967,0,0.7745967 };
        const double Ci[3] = { 0.5555556,0.8888889,0.5555556 };


        double Q, S = 0.0;
        for (int i = 0; i < 3; i++)
        {
            Q = (p1 + p2) / 2 + ((p2 - p1) / 2) * Xi[i];
            S += Ci[i] * A.GetRes(Q);
        }
        return (p2 - p1) / 2 * S;
    }
    double SquarePm(double i)
    {
        return(A.GetRes(i) * h);
    }
    double SquarePpogr(double i)
    {
        return(2 * abs(powf(A.GetRes(i + h) - A.GetRes(i),3) / (2*h)));
    }
    double SquareTm(double i)
    {
        return (h * ((A.GetRes(i - h) + A.GetRes(i))) / 2);
    }
    double SquareTpogr(double i, double max)
    {
        if (A.GetRes(i) > max)
        {
            max = A.GetRes(i);
        }
        return(max);
    }
    double Parabm(double i, double h)
    {
        return((h / 6) * (A.GetRes(i - h) + A.GetRes(i) + 4 * A.GetRes((2 * i - h) / 2)));
    }
    double Montem(double i,double h)
    {
        double a = (i / h + (rand() % (int)(1/h))*h) * h;
        return (A.GetRes(a));
    }
    double Runge6(double i, double h)
    {
        double k1, k2, k34, k5, k6;
        double a = 0;
        k1 = A.GetRes(i);
        k2 = A.GetRes(i + h / 4);
        k34 = A.GetRes(i + h / 2);
        k5 = A.GetRes(i + 3 * h / 4);
        k6 = A.GetRes(i + h);
        a = (7.0 / 90 * (k1 + k6) + 16.0 / 45 * (k2 + k5) + 2.0 / 15 * k34) * h;
        return a;
    }
    double Runge3(double i, double h)
    {
        double k1;
        double k2;
        double a = 0;
        k1 = A.GetRes(i);
        k2 = A.GetRes(i + 2 * h / 3);
        a = (h / 4) * (k1 + 3 * k2);
        return a;
    }
    void Round()
    {
        Res1.Integ = round(Res1.Integ*(1/toch))*toch;
        Res1.pogr = round(Res1.pogr*(1/toch))*toch;
    }
public:
    Integral(double x1, double x2, std::string str_)
    {
        h = 0.01;
        p1 = x1;
        p2 = x2;
        str = str_;
        A.SetStr(str);
        toch = 0.001;
    }
    Integral(double x1, double x2, double h_, std::string str_)
    {
        h = h_;
        p1 = x1;
        p2 = x2;
        str = str_;
        A.SetStr(str);
        toch = 0.001;

    }
    void set_str(std::string a)
    {
        str = a;
    }
    void set_toch(double a)
    {
        toch = a;
    }
    Res SquareP()
    {

        Res1.Integ = 0;
        Res1.pogr = 0;
        for (double i = p1; i <= p2; i = i + h)
        {
            Res1.Integ += SquarePm(i);
            Res1.pogr += SquarePpogr(i);
        }
        Round();
        return(Res1);
    }
    Res SquareT()
    {
        Res1.Integ = 0;
        double max = 0;
        for (double i = p1 + h; i <= p2; i = i + h)
        {
            Res1.Integ += SquareTm(i);
            max = SquareTpogr(i, max);
            Res1.pogr = (p2 - p1) * max * h / 2;
        }
        Round();
        return(Res1);
    }
    Res Parab()
    {
        Res1.Integ = 0;
        double a = 0;
        for (double i = p1 + h; i <= p2; i = i + h)
        {
            Res1.Integ += Parabm(i, h);
        }
        for (double i = p1 + h / 2; i <= p2; i = i + h / 2)
        {
            a += Parabm(i, h / 2);
        }
        Res1.pogr = 2 * abs(a - Res1.Integ);
        Round();
        return(Res1);
    }
    Res Gauss()
    {
        Res1.Integ = 0;
        double a = 0;
        for (size_t i = 0; i < (p2 - p1) / h; ++i)
        {
            Res1.Integ += Gaussm(p1 + i * h, p1 + (i + 1) * h);
        }
        for (size_t i = 0; i < (p2 - p1) / (h / 2); ++i)
        {
            a += Gaussm(p1 + i * h / 2, p1 + (i + 1) * h / 2);
        }
        Res1.pogr = abs(a - Res1.Integ) * 2;
        Round();
        return(Res1);
    }
    Res Monte()
    {
        Res1.Integ = 0;
        srand(time(0));
        double a = 0;
        for (double i = p1; i < p2; i=i+h)
        {
            Res1.Integ += Montem(i,h);
        }
        for (double i = p1; i < p2; i=i+h/2)
        {
            a += Montem(i,h/2);
        }
        Res1.Integ = Res1.Integ * ((double)(p2 - p1) / (int)((p2 - p1) / h));
        a = a * ((double)(p2 - p1) / (int)(2 * ((p2 - p1) / h)));
        Res1.pogr = 2 * abs(Res1.Integ - a);
        Round();
        return(Res1);
    }
    Res Runge(int n)
    {
        Res1.Integ = 0;
        Res1.pogr = 0;
        double a = 0;
        if (n == 2)
        {
            Res1 = SquareT();
        }
        if (n == 3)
        {
            for (double i = p1; i <= p2; i = i + h)
            {
                Res1.Integ += Runge3(i, h);
            }
            for (double i = p1; i <= p2; i = i + h / 2)
            {
                a += Runge3(i, h / 2);
            }
            Res1.pogr = 2 * abs(Res1.Integ - a);
            Round();
        }
        if (n == 4)
        {
            Res1 = Parab();

        }
        if (n == 6)
        {
            for (double i = p1; i <= p2; i = i + h)
            {
                Res1.Integ += Runge6(i, h);
            }
            for (double i = p1; i <= p2; i = i + h / 2)
            {
                a += Runge6(i, h / 2);
            }
            Res1.pogr = 2 * abs(Res1.Integ - a);
            Round();
        }
 

  
        
        return (Res1);
    }
};
std::ostream& operator<< (std::ostream& out, Res Res1)
{
    out << Res1.Integ << "+/-" << Res1.pogr;
    return (out);
}
int main(int argc, char** argv)
{
    const char* input = "sinx";
    double x1 = 0;
    double x2 = 3.14;
    double step = 0.01;
    if (argc > 1)
        input = argv[1];
    if (argc > 2)
        x1 = strtod(argv[2], nullptr);
    if (argc > 3)
        x2 = strtod(argv[3], nullptr);
    if (argc > 4)
        step = strtod(argv[4], nullptr);

    double p;
    Integral<double> a(x1, x2, step, input);
    std::cout << a.SquareP() << std::endl;
    std::cout << a.SquareT() << std::endl;
    std::cout << a.Parab() << std::endl;
    std::cout << a.Gauss() << std::endl;
    std::cout << a.Monte() << std::endl;
    std::cout << a.Runge(6) << std::endl;
    return 0;
}