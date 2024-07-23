#pragma once
#include "matrix.h"
#include "substance.h"
#include <list>
#include <map>
#include <regex>
#include <sstream>
#include <numeric>
#include <string>
using namespace std;

struct condition {
	__readonly double t;//开尔文温度
	__readonly double p;//反应压强

	constexpr condition(double t_, double p_) :t(t_), p(p_) {
		if (t_ < 0 || p_ < 0)throw L"无效的温度或压强";
	}
	constexpr condition(condition&&)noexcept = default;
	constexpr condition(const condition&)noexcept = default;
};
const constexpr condition stp(273.15 + 25, 101.00);//0C->K 1atm->kPa

class equation {
private:
	list<pair<substance, bool>>substances;//物质列表
	condition tp;//温度&压强
	ratmatrix mat;//储存系数结果矩阵 外层列 元素|电荷 内层行 未知数|系数

public:
	equation(const wstring& str, const condition& cond) :tp(cond) {
		//按连接符号分割 判断空
		size_t p = str.find(L"--");
		if (p == wstring::npos || p != str.rfind(L"--"))
			throw wstring(L"连接符号错误");
		if (p == 0 || p + 2 == str.size())throw wstring(L"化学式不能为空");
		//分为左右，正则进一步分割
		const wstring::const_iterator split_begin = str.cbegin() + p, split_end = split_begin + 2;
		static const wregex split_reg(L"\\+(?!\\>)");
		static const wsregex_token_iterator it_end;
		size_t pos = 0;
		try {
			for (wsregex_token_iterator it(str.cbegin(), split_begin, split_reg, -1); it != it_end; it++) {
				pos++;
				substances.emplace_back(substance(it->first, it->second), false);
			}
			for (wsregex_token_iterator it(split_end, str.cend(), split_reg, -1); it != it_end; it++) {
				pos++;
				substances.emplace_back(substance(it->first, it->second), true);
			}
		}
		catch (const wchar_t* pt) {
			throw L"物质" + to_wstring(pos) + L':' + pt;
		}
	}

	void balance(bool int_ = true) {//是否整数化
		//统计元素及其个数
		map<short, rat>elems;//统计每种元素
		for (const auto& pair_ : substances)
			for (const auto& p : pair_.first.elements)
				elems[p.first] += pair_.first.count * p.second * (pair_.second ? -1 : 1);
		//判断有无待定系数
		bool has = false;
		for (const auto& p : substances)
			if (p.first.count != 0) {
				has = true;
				break;
			}
		//向矩阵中插入元素
		for (const auto& pair_ : elems) {
			vector<rat> nums;
			for (const auto& p : substances)
				if (p.first.count == 0) {
					const auto it = p.first.elements.find(pair_.first);
					nums.push_back(it != p.first.elements.cend() ? it->second * (p.second ? -1 : 1) : 0);
				}
			if (has)nums.push_back(pair_.second);
			mat.append(nums);
		}
		//求解赋值
		if (mat.sizev() == 1)throw wstring(L"所有系数已被决定");
		mat.to_upper_triangular();
		if (!mat.solve())throw wstring(L"矩阵无解");
		//把待定系数项放回去
		if (has) {
			ratmatrix mat_;
			int i = 0;
			for (const auto& p : substances) {
				if (p.first.count == 0) {
					mat_.append(mat[i]);
					i++;
					continue;
				}
				vector<rat>vec = mat.back();
				for (rat& r : vec)r *= p.first.count;
				mat_.append(move(vec));
			}
			mat = move(mat_);
		}
		//整数化
		if (int_) {
			for (size_t i = 0; i < mat.sizev(); i++) {
				int x = 0, y = 1;//分子，分母
				for (int j = 0; j < mat.sizeh(); j++) {
					x = int(gcd(mat[j][i].up, x));
					y = int(lcm(mat[j][i].down, y));
				}
				rat r(x, y);
				for (size_t j = 0; j < mat.sizeh(); j++)
					mat.in(j, i) /= r;
			}
		}
	}

	[[nodiscard]] wstring print()const {//输出配平结果
		auto join = [&](const list<wstring>& l, const wchar_t* w) {
			if (l.size() == 1)return l.front();
			int i = 0;
			wstring result;
			for (const wstring& str : l) {
				result.append(str);
				result.append(w);
				i++;
				if (i == l.size() - 1)break;
			}
			result.append(l.back());
			return result;
			};
		auto to_wstring_ = [](double num, unsigned base) {
			wostringstream out;
			out.precision(base);
			out << fixed << num;
			return out.str();
			};
		wstring cond = tp.t == stp.t ? L"" : to_wstring_(tp.t, 2) + L"K";
		if (tp.p != stp.p)cond += (cond.empty() ? L"" : L",") + to_wstring_(tp.p, 2) + L"kPa";
		cond = L"==" + cond + L"==";
		list<wstring> list_;
		for (int i = 0; i < mat.sizev(); i++) {//遍历列
			list<wstring> left, right;
			int j = -1;
			for (const pair<substance, bool>& p : substances) {
				j++;
				if (mat[j][i] == 0)continue;
				wstring text = (abs(mat[j][i]) != 1 ? L"<span style=\"background-color: #3c9eef;\">"
					+ (abs(mat[j][i])).to_wstring() + L"</span>" : L"") + p.first.html;
				if (mat[j][i] > 0 != p.second) { left.push_back(move(text)); }
				else { right.push_back(move(text)); }
			}
			list_.push_back(join(left, L" + ") + cond + join(right, L" + "));
		}
		list<wstring> empty;
		int i = -1;
		for (const pair<substance, bool>& p : substances) {
			i++;
			if (mat.count(i, 0) == mat.sizev()) { empty.push_back(p.first.html); }
		}
		if (!empty.empty()) { list_.push_back(L"可能未参与物质:" + join(empty, L" ")); }
		return join(list_, L"<br>");
	}

	inline equation(const equation&)noexcept = default;
	inline equation(equation&&)noexcept = default;
	inline ~equation() = default;
};