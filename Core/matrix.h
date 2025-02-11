#include "rat.h"
#include <cassert>
#include <iosfwd>
#include <set>
#include <vector>

template<std::common_with<short> _Num>
class basic_matrix {//未说明默认行操作
private:
	std::vector<std::vector<_Num>>mat;

public:
	using const_iterator = std::vector<std::vector<_Num>>::const_iterator;//行只读迭代器
	static constexpr auto npos{ static_cast<size_t>(-1) };

	constexpr basic_matrix() = default;
	constexpr basic_matrix(const basic_matrix&) = default;
	constexpr basic_matrix(basic_matrix&&) = default;
	constexpr ~basic_matrix() = default;
	constexpr void operator =(const basic_matrix& x) { mat = x.mat; }
	constexpr void operator =(basic_matrix&& x)noexcept { mat = x.mat; }
	[[nodiscard]] constexpr auto operator<=>(const basic_matrix& x)const { return mat <=> x.mat; }
	[[nodiscard]] constexpr bool operator==(const basic_matrix& x)const { return x.mat == mat; }

	[[nodiscard]] constexpr size_t sizeh()const { return mat.size(); }//行数
	[[nodiscard]] constexpr size_t sizev()const { return mat.empty() ? 0 : mat.front().size(); }//列数
	[[nodiscard]] constexpr bool empty()const { return mat.empty(); }
	[[nodiscard]] constexpr const std::vector<_Num>& operator[](const size_t& x)const { return mat[x]; }
	[[nodiscard]] constexpr _Num& in(const size_t& x, const size_t& y) { return mat[x][y]; }
	[[nodiscard]] constexpr _Num& in(const_iterator x, const size_t& y) { return *x[y]; }
	[[nodiscard]] constexpr const std::vector<_Num>& front()const { return mat.front(); }
	[[nodiscard]] constexpr const std::vector<_Num>& back()const { return mat.back(); }
	[[nodiscard]] constexpr const_iterator cbegin()const { return mat.cbegin(); }
	[[nodiscard]] constexpr const_iterator cend()const { return mat.cend(); }
	[[nodiscard]] constexpr const_iterator begin()const { return mat.cbegin(); }
	[[nodiscard]] constexpr const_iterator end()const { return mat.cend(); }

	constexpr void append(const std::vector<_Num>& v) {
		assert(!(v.empty() || !mat.empty() && v.size() != sizev()));
		mat.push_back(v);
	}

	constexpr void append(std::vector<_Num>&& v) {//通过move元素来构造
		assert(!(v.empty() || !mat.empty() && v.size() != sizev()));
		mat.push_back(v);
	}

	constexpr void transpose()noexcept {//行列交换
		if (mat.empty())return;
		std::vector<std::vector<_Num>>_mat(sizev(), vector<_Num>(sizeh()));
		for (size_t i = 0; i < sizeh(); i++)
			for (size_t j = 0; j < sizev(); j++)
				_mat[j][i] = std::move(mat[i][j]);
		mat = std::move(_mat);
	}

	template<typename _Char, typename _Traits = std::char_traits<_Char>>
		requires std::_Is_character<_Char>::value
	friend inline std::basic_ostream<_Char, _Traits>& operator<<//debugging support
		(std::basic_ostream<_Char, _Traits>& os, const basic_matrix& m) {
		os << os.widen('{');
		for (const auto& v : m.mat) {
			os << os.widen('{');
			for (const auto& num : v) {
				os << num << os.widen(',');
			}
			os << os.widen('}');
		}
		os << os.widen('}');
		return os;
	}

	constexpr void operator *=(const _Num& t)noexcept {
		for (auto& v : mat) {
			for (_Num& _t : v) {
				_t *= t;
			}
		}
	}

	inline void to_upper_triangular() {//化为上三角阵
		if (mat.empty())return;
		for (size_t i = 0; i < sizev(); i++) {//i 列
			for (size_t j = i + 1; j < sizeh(); j++) {//j行
				if (mat.at(j).at(i) == 0)continue;
				std::swap(mat[i], mat[j]);
				if (mat.at(j).at(i) == 0)continue;
				_Num a = mat.at(j).at(i) / mat.at(i).at(i);
				for (size_t k = i; k < sizev(); k++) {//倍增后减过去
					if (mat.at(i).at(k) != 0) {
						mat[j][k] -= mat.at(i).at(k) * a;
					}
				}
			}
		}
		while (sizeh() > sizev())mat.pop_back();//删除空白末尾
	}

	[[nodiscard]] bool solve() {
		if (mat.empty())return false;
		std::set<size_t>pivot_columns;
		//化为简化行阶梯阵，移除空行，获得主元列
		for (size_t i = sizeh() - 1; i != npos; i--) {
			size_t pivot_col = i;
			while (pivot_col < sizev() && mat.at(i).at(pivot_col) == 0)
				pivot_col++;
			if (pivot_col == sizev()) {
				mat.erase(mat.begin() + i);
				continue;
			}
			pivot_columns.insert(pivot_col);//获得主元列
			//系数化1
			if (mat.at(i).at(pivot_col) != 1) {
				_Num factor = mat.at(i).at(pivot_col);
				for (size_t j = pivot_col; j < sizev(); j++) {
					mat[i][j] /= factor;
				}
			}
			//消除主元列其他元素
			for (size_t j = 0; j < i; j++) {
				if (mat.at(j).at(pivot_col) == 0)continue;
				_Num factor = mat.at(j).at(pivot_col);
				for (size_t k = pivot_col; k < sizev(); k++) {
					mat[j][k] -= mat.at(i).at(k) * factor;
				}
			}
		}
		//判断无解情况
		if (sizeh() >= sizev() - 1)
			for (const auto& line : mat)
				if (line.cend() - 1 == find_if_not(line.cbegin(), line.cend(), [](const rat& x) {return x == 0; }))
					return false;
		//生成-F矩阵 拼接I矩阵
		size_t original_size = sizev();//原本列数
		for (size_t i = sizev() - 1; i != npos; i--)
			if (pivot_columns.contains(i))
				for (size_t j = 0; j < sizeh(); j++)
					mat[j].erase(mat[j].begin() + i);
		*this *= -1;
		for (size_t i = 0; i < sizev(); i++) {
			std::vector<_Num>vec(sizev(), 0);
			vec[i] = 1;
			mat.push_back(move(vec));
		}
		//重新排序
		if (sizev() == 1)return true;
		std::vector<size_t>order;
		for (const size_t& x : pivot_columns)
			order.push_back(x);
		for (size_t i = 0; i < original_size; i++)
			if (!pivot_columns.contains(i))order.push_back(i);
		std::vector<std::vector<_Num>>reordered_mat;
		for (size_t i = 0; i < original_size; i++)
			reordered_mat.push_back(move(mat[order.at(i)]));
		mat = move(reordered_mat);
		return true;
	}
};

using ratmatrix = basic_matrix<rat>;