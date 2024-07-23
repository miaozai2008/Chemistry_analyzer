#pragma once
#include <cassert>
#include <string>
#include <type_traits>
#include <iosfwd>
#include <numeric>

template<typename _Elem>
class basic_rat {
	static_assert(std::is_signed<_Elem>::value&& std::is_integral<_Elem>::value, "_Elem must be a signed integer");

private:
	constexpr void simplify() {
		assert(down);
		if (up == 0) {
			down = 1;
			return;
		}
		const _Elem x = _Elem(std::gcd(up, down));
		up /= x;
		down /= x;
		if (down < 0) {
			up = -up;
			down = -down;
		}
	}

public:
	__readonly _Elem up;//分子
	__readonly _Elem down;//分母

	constexpr basic_rat(const _Elem& x, const _Elem& y) :up(x), down(y) { simplify(); }
	constexpr basic_rat(_Elem&& x, _Elem&& y) : up(x), down(y) { simplify(); }
	constexpr basic_rat(const _Elem& x = 0) : up(x), down(1) {  }
	constexpr basic_rat(const basic_rat&)noexcept = default;
	constexpr basic_rat(basic_rat&&)noexcept = default;
	constexpr ~basic_rat() = default;

	constexpr void operator=(const basic_rat& x)noexcept { up = x.up; down = x.down; }
	constexpr void operator=(basic_rat&& x)noexcept { up = x.up; down = x.down; }
	[[nodiscard]] constexpr basic_rat operator-()const noexcept { return basic_rat(-up, down); }
	friend std::ostream& operator<<(std::ostream& os, const basic_rat& x) { os << x.to_string(); return os; }
	friend std::wostream& operator<<(std::wostream& os, const basic_rat& x) { os << x.to_wstring(); return os; }

	[[nodiscard]] constexpr friend basic_rat abs(const basic_rat& r)noexcept { return r.up < 0 ? -r : r; }
	[[nodiscard]] constexpr friend basic_rat abs(basic_rat&& r)noexcept { return r.up < 0 ? -r : r; }
	[[nodiscard]] constexpr double to_double()const { return double(up) / down; }
	[[nodiscard]] constexpr std::string to_string()const { return std::to_string(up) + (down == 1 ? "" : "/" + std::to_string(down)); }
	[[nodiscard]] constexpr std::wstring to_wstring()const { return std::to_wstring(up) + (down == 1 ? L"" : L"/" + std::to_wstring(down)); }

	[[nodiscard]] constexpr bool operator==(const basic_rat& x)const noexcept { return down == x.down && up == x.up; }
	[[nodiscard]] constexpr bool operator!=(const basic_rat& x)const noexcept { return down != x.down || up != x.up; }
	[[nodiscard]] constexpr bool operator>(const basic_rat& x)const noexcept { return up * x.down > x.up * down; }
	[[nodiscard]] constexpr bool operator<(const basic_rat& x)const noexcept { return up * x.down < x.up * down; }
	[[nodiscard]] constexpr bool operator>=(const basic_rat& x)const noexcept { return *this == x || *this > x; }
	[[nodiscard]] constexpr bool operator<=(const basic_rat& x)const noexcept { return *this == x || *this < x; }

	[[nodiscard]] constexpr basic_rat operator*(const basic_rat& x)const { return basic_rat(x.up * up, x.down * down); }
	constexpr void operator*=(const basic_rat& x) { *this = basic_rat(x.up * up, x.down * down); }
	[[nodiscard]] constexpr basic_rat operator/(const basic_rat& x)const { return basic_rat(x.down * up, x.up * down); }
	constexpr void operator/=(const basic_rat& x) { *this = basic_rat(x.down * up, x.up * down); }
	[[nodiscard]] constexpr basic_rat operator+(const basic_rat& x)const { return basic_rat(up * x.down + x.up * down, down * x.down); }
	constexpr void operator+=(const basic_rat& x) { *this = basic_rat(up * x.down + x.up * down, down * x.down); }
	[[nodiscard]] constexpr basic_rat operator-(const basic_rat& x)const { return basic_rat(up * x.down - x.up * down, down * x.down); }
	constexpr void operator-=(const basic_rat& x) { *this = basic_rat(up * x.down - x.up * down, down * x.down); }
};

using rat = basic_rat<int>;