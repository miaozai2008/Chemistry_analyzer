#pragma once
#include "matrix.h"
#include "substance.h"
#include <list>
#include <map>
#include <numeric>
#include <regex>
#include <sstream>
#include <string>
using namespace std;

class equation {
public:
	struct condition {
		double t;//开尔文温度
		double p;//反应压强

		constexpr condition(condition&&)noexcept = default;
		constexpr condition(const condition&)noexcept = default;
		constexpr condition(double t_, double p_) :t(t_), p(p_) {
			if (t_ < 0 || p_ < 0)throw L"无效的温度或压强";
		}
	};

private:
	list<pair<substance, bool>>substances;//物质列表
	condition tp;//温度&压强
	ratmatrix mat;//储存系数结果矩阵 外层列 元素|电荷 内层行 未知数|系数
	bool has_more = true;//是否可以计算熵焓

public:
	equation(const wstring& str, condition&& cond) :tp(cond) {
		//按连接符号分割 判断空
		size_t p = str.find(L"--");
		if (p == wstring::npos || p != str.rfind(L"--"))throw wstring(L"连接符号错误");
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

	void compute(bool int_ = true) {//是否整数化 配平 计算h s g
		//统计元素及其个数
		map<short, rat>elems;//统计每种元素
		for (const auto& [sub, reversed] : substances)
			for (const auto& [elem, num] : sub.elements)
				elems[elem] += sub.count * num * (reversed ? -1 : 1);
		//判断有无待定系数
		bool has = false;
		for (const auto& [sub, _] : substances)
			if (sub.count != 0) {
				has = true;
				break;
			}
		//向矩阵中插入元素
		for (const auto& [elem, num] : elems) {
			vector<rat>nums;
			for (const auto& [sub, reversed] : substances)
				if (sub.count == 0) {
					const auto it = sub.elements.find(elem);
					nums.push_back((it == sub.elements.cend() ? 0 : it->second) * (reversed ? -1 : 1));
				}
			if (has)nums.push_back(num);
			mat.append(move(nums));
		}
		//求解赋值
		if (mat.sizev() == 1)throw wstring(L"所有系数已被决定");
		mat.to_upper_triangular();
		if (!mat.solve())throw wstring(L"矩阵无解");
		//把待定系数项放回去
		if (has) {
			ratmatrix mat_;
			size_t i = 0;
			for (const auto& [sub, _] : substances) {
				if (sub.count == 0) {
					mat_.append(mat[i]);
					i++;
					continue;
				}
				vector<rat>vec = mat.back();
				for (rat& r : vec)r *= sub.count;
				mat_.append(move(vec));
			}
			mat = move(mat_);
		}
		//整数化
		if (int_) {
			for (size_t i = 0; i < mat.sizev(); i++) {
				int x = 0, y = 1;//分子，分母
				for (int j = 0; j < mat.sizeh(); j++) {
					x = gcd(mat[j][i].up, x);
					y = lcm(mat[j][i].down, y);
				}
				rat r(x, y);
				for (size_t j = 0; j < mat.sizeh(); j++)
					mat.in(j, i) /= r;
			}
		}
		//计算s h g
		for (auto& [sub, _] : substances) {
			if (!sub.search()) {
				has_more = false;
				return;
			}
		}
	}

	[[nodiscard]] wstring print()const noexcept {//输出结果
		auto join = [&](const list<wstring>& l, const wchar_t* w) {
			if (l.empty())return wstring();
			if (l.size() == 1)return l.front();
			size_t i = 0;
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
		wstring cond = tp.t == 273.15 + 25 ? L"" : to_wstring_(tp.t, 2) + L"K";//stp t
		if (tp.p != 101)cond += (cond.empty() ? L"" : L",") + to_wstring_(tp.p, 2) + L"kPa";//stp p
		cond = L"==" + cond + L"==";
		list<wstring> list_;
		for (size_t i = 0; i < mat.sizev(); i++) {//遍历列
			list<wstring> left, right;
			int j = -1;
			for (const auto& [sub, reversed] : substances) {
				j++;
				if (mat[j][i] == 0)continue;
				wstring text = (abs(mat[j][i]) != 1 ? L"<span style=\"background-color: #3c9eef;\">"
					+ (abs(mat[j][i])).to_wstring() + L"</span>" : L"") + sub.html;
				if (mat[j][i] > 0 != reversed)  left.push_back(move(text));
				else right.push_back(move(text));
			}
			list_.push_back(join(left, L" + ") + cond + join(right, L" + "));
			j = -1;
			if (has_more) {//添加S H G
				double h = 0, s = 0;
				for (const auto& [sub, reversed] : substances) {
					j++;
					h += sub.h * mat[j][i].to_double() * (reversed ? 1 : -1);
					s += sub.s * mat[j][i].to_double() * (reversed ? 1 : -1);
				}
				list_.push_back(L" ΔH=<u>" + to_wstring_(h, 3) + L"</u>kJ/mol");
				if (tp.p == 101) {//标准大气压
					list_.back().append(L" ΔS=<u>" + to_wstring_(s, 3) + L"</u>J/(molK)");
					list_.back().append(L" ΔG=<u>" + to_wstring_(h - s * tp.t / 1000, 3) + L"</u>kJ/mol");
				}
			}
		}
		list<wstring> empty;
		int i = -1;
		for (const auto& [sub, _] : substances) {
			i++;
			if (mat.count(i, 0) == mat.sizev())empty.push_back(sub.html);
		}
		if (!empty.empty())list_.push_back(L"可能未参与物质:" + join(empty, L" "));
		if (!has_more)list_.push_back(L"无算计算熵焓:未指定物态或未搜寻到数据");
		return join(list_, L"<br>");
	}

	inline equation(const equation&)noexcept = default;
	inline equation(equation&&)noexcept = default;
	inline ~equation() = default;
};