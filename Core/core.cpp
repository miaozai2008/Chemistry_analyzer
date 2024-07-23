#include "equation.h"
#include <iostream>
#include <locale>
#include <string>
using namespace std;

int main() {
	wcout.imbue(locale("zh_CN", LC_CTYPE));//处理wcout的错误
	wcin.imbue(locale("zh_CN", LC_CTYPE));//处理wcin错误
	wstring formula;
	double tem = 0, pres = 0;
	while (true) {
		wcin >> formula;
		if (formula == L"[exit]")break;
		wcin >> tem >> pres;
		try {
			equation my_equation(formula, condition(tem, pres));
			my_equation.balance();
			wcout << my_equation.print() << endl;
		}
		catch (const wstring& error) {
			wcout << L"[错误]" << error << endl;
		}
	}
	return 0;
}