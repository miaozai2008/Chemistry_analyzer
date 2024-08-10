#include "rat.h"
#include <fstream>
#include <map>
#include <regex>
#include <stack>
#include <string>
#include <unordered_map>
using namespace std;

class substance {
private:
	wstring formula;//物质化学式
	using _Cwit = wstring::const_iterator;

	[[nodiscard]] inline map<short, int>dispose(_Cwit begin, _Cwit end)const {
		if (!iswupper(*begin))throw L"元素必须以大写开头";
		map<short, int>elems;
		stack<int>amps;
		int amp_ = 1;
		wstring num;
		for (_Cwit it = end - 1;; it--) {
			if (iswdigit(*it))num.append(&*it);
			else if (*it == L'(' || *it == L'[') {
				amp_ /= amps.top();
				amps.pop();
			}
			else {
				const int cnt = num.empty() ? 1 : stoi(num);
				num.clear();
				amps.push(cnt);
				amp_ *= cnt;
				if (*it == L')' || *it == L']')continue;
				if (iswupper(*it))elems[*it - L'A'] += amp_;
				else if (iswlower(*it)) {
					short key = (*it - L'a' + 1) * 26;
					it--;
					if (!iswupper(*it))throw L"无效元素";
					elems[*it - L'A' + key] += amp_;
				}
				else throw L"未知符号";
				amp_ /= cnt;
				amps.pop();
			}
			if (it == begin)return elems;
		}
	}

public:
	__readonly wstring html;//对外显示化学式
	__readonly rat count = 0;//用户定义系数
	__readonly map<short, int> elements; //元素及个数列表
	__readonly double h = 0;//焓 kJ
	__readonly double s = 0;//熵 J

	substance(_Cwit begin, _Cwit end) {
		wsmatch match;
		//分离状态
		wstring state;
		static const wregex form_reg(L"(\\([a-z]+\\))$");
		bool has_state = true;
		if (regex_search(begin, end, match, form_reg)) {
			state = match[1];
			end -= match.length();
		}
		else has_state = false;
		//是否为空 检查括号是否匹配 数字超过8位
		if (begin == end)throw L"物质不能为空";
		else {
			int x = 0, y = 0, cnt = 0;
			for (_Cwit it = begin; it < end; it++)
				if (iswdigit(*it)) {
					if (*it == L'0' && cnt == 0)throw L"非法数字";
					else cnt++;
				}
				else {
					if (*it == L'(')x++;
					else if (*it == L')')x--;
					else if (*it == L'[') {
						if (x > 0)throw L"括号嵌套错误";
						y++;
					}
					else if (*it == L']') {
						if (x > 0)throw L"括号嵌套错误";
						y--;
					}
					if (cnt > 8 || x < 0 || y < 0)throw L"数字过大或括号不匹配";
					cnt = 0;
				}
			if (x != 0 || y != 0 || cnt > 8)throw L"数字过大或括号不匹配";
		}
		//分离系数
		static const wregex rat_reg(L"^(\\d+)/?(\\d+)?");
		if (regex_search(begin, end, match, rat_reg)) {
			count = rat(stoi(match[1]), (match.length(2) ? stoi(match[2]) : 1));
			begin += match.length();
		}
		//分离电荷数
		int cnum = 0;
		static const wregex charge_reg(L"<(\\d*)(\\+|\\-)>$");
		if (regex_search(begin, end, match, charge_reg)) {
			cnum = (match.length(1) ? stoi(match[1]) : 1) * (match[2] == L'-' ? -1 : 1);
			end -= match.length();
		}
		//填写formula
		if (has_state)
			formula = wstring(begin, end) + (cnum ? (cnum > 0 ? L"+" : L"") + to_wstring(cnum) : L"") + state;
		//判断是否空 ·的存在 临时变量
		if (begin == end)throw L"物质不能为空";
		html.assign(begin, end);
		bool has_dot = html.find(L'·') != wstring::npos;
		//分离出元素和对应个数 支持·
		if (!has_dot)elements = dispose(begin, end);
		else {
			if (*begin != L'[' || *(end - 1) != L']')throw L"物质格式错误，首尾应有[]";
			begin++;
			end--;
			static const wregex split_reg(L"·"), num_reg(L"^(\\d*)");
			static const wsregex_token_iterator it_end;
			for (wsregex_token_iterator it(begin, end, split_reg, -1); it != it_end; it++) {
				regex_search(it->first, it->second, match, num_reg);
				const int num = match[1].length() ? stoi(match[1]) : 1;
				if (match[1].second == it->second)throw L"物质格式错误，不能为空";
				for (const auto& p : dispose(match[1].second, it->second))
					elements[p.first] += num * p.second;
			}
		}
		if (cnum != 0)elements.insert(pair<short, int>(-1, cnum));
		//将表达式html化
		for (size_t i = 1; i < html.size(); i++) {
			if (iswdigit(html.at(i)) && html.at(i - 1) != L'·') {
				html.insert(i, L"<sub>");
				i += 5;
				for (; i < html.size(); i++)if (!iswdigit(html.at(i)))break;
				html.insert(i, L"</sub>");
				i += 6;
			}
		}
		if (cnum != 0) {
			if (abs(cnum) == 1)  html += L"<sup>" + wstring(cnum < 0 ? L"-" : L"+") + L"</sup>";
			else  html += L"<sup>" + to_wstring(abs(cnum)) + (cnum < 0 ? L"-" : L"+") + L"</sup>";
		}
		if (has_dot)html = L'[' + html + L']';
		if (has_state)html += state;
	}

	[[nodiscard]] bool search() {//查询txt文件找熵焓
		if (formula.empty())return false;
		static unordered_map<wstring, pair<double, double>> datas;//物质 H S
		if (datas.empty()) {
			wifstream input("datas.txt");
			if (!input.is_open())throw wstring(L"datas.txt未找到");
			wstring line;
			wregex reg(L"(.+?)\t(.+?)\t(.+?)");
			wsmatch match;
			while (getline(input, line)) {
				if (!regex_match(line, match, reg))throw wstring(L"datas.txt格式错误");
				datas.emplace(match[1], make_pair(stod(match[2]), stod(match[3])));
			}
		}
		decltype(datas)::const_iterator it = datas.find(formula);
		if (it == datas.cend())return false;
		h = it->second.first;
		s = it->second.second;
		return true;
	}

	inline substance(const substance&)noexcept = default;
	inline substance(substance&&)noexcept = default;
	inline ~substance() = default;
};