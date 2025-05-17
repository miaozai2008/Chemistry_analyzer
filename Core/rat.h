#pragma once
#include <concepts>
#include <numeric>
#include <stdexcept>
#include <string>
#include <compare>

template<std::signed_integral _Num>
class basic_rat {
private:
	constexpr void simplify() {
		if (!down)throw std::logic_error("rat: divide by 0");
		if (up == 0) {
			down = 1;
			return;
		}
		const _Num x = std::gcd(up, down);
		up /= x;
		down /= x;
		if (down < 0) {
			up = -up;
			down = -down;
		}
	}

public:
	__readonly _Num up;//分子
	__readonly _Num down;//分母

	constexpr basic_rat(const _Num& x, const _Num& y) : up(x), down(y) { simplify(); }
	constexpr basic_rat(const _Num& x = 0) : up(x), down(1) {}
	constexpr basic_rat(const basic_rat&) noexcept = default;
	constexpr basic_rat(basic_rat&&) noexcept = default;
	constexpr ~basic_rat() = default;

	constexpr basic_rat& operator=(const basic_rat& x) noexcept {
		if (this != &x) {
			up = x.up;
			down = x.down;
		}
		return *this;
	}
	constexpr basic_rat& operator=(basic_rat&& x) noexcept {
		if (this != &x) {
			up = x.up;
			down = x.down;
		}
		return *this;
	}

	[[nodiscard]] constexpr basic_rat operator-()const noexcept { return basic_rat(-up, down); }
	friend std::ostream& operator<<(std::ostream& os, const basic_rat& x) { os << x.to_string(); return os; }
	friend std::wostream& operator<<(std::wostream& os, const basic_rat& x) { os << x.to_wstring(); return os; }
	[[nodiscard]] constexpr std::strong_ordering operator<=>(const basic_rat& x)const noexcept { return up * x.down <=> x.up * down; }
	[[nodiscard]] constexpr std::strong_ordering operator<=>(const _Num& x)const noexcept { return up <=> x * down; }
	[[nodiscard]] constexpr bool operator==(const _Num& x)const noexcept { return down == 1 && up == x; }
	[[nodiscard]] constexpr bool operator==(const basic_rat& x)const noexcept { return up == x.up && down = x.down; }

	[[nodiscard]] constexpr friend basic_rat abs(const basic_rat& r)noexcept { return r.up < 0 ? -r : r; }
	[[nodiscard]] constexpr double to_double()const { return double(up) / down; }
	[[nodiscard]] constexpr std::string to_string()const { return std::to_string(up) + (down == 1 ? "" : "/" + std::to_string(down)); }
	[[nodiscard]] constexpr std::wstring to_wstring()const { return std::to_wstring(up) + (down == 1 ? L"" : L"/" + std::to_wstring(down)); }

	[[nodiscard]] constexpr basic_rat operator*(const basic_rat& x)const noexcept { return { x.up * up, x.down * down }; }
	constexpr basic_rat& operator*=(const basic_rat& x)noexcept { return *this = { x.up * up, x.down * down }; }
	[[nodiscard]] constexpr basic_rat operator/(const basic_rat& x)const { return { x.down * up, x.up * down }; }
	constexpr basic_rat& operator/=(const basic_rat& x) { return *this = { x.down * up, x.up * down }; }
	[[nodiscard]] constexpr basic_rat operator+(const basic_rat& x)const noexcept { return { up * x.down + x.up * down, down * x.down }; }
	constexpr basic_rat& operator+=(const basic_rat& x) noexcept { return *this = { up * x.down + x.up * down, down * x.down }; }
	[[nodiscard]] constexpr basic_rat operator-(const basic_rat& x)const noexcept { return { up * x.down - x.up * down, down * x.down }; }
	constexpr basic_rat& operator-=(const basic_rat& x)noexcept { return *this = { up * x.down - x.up * down, down * x.down }; }
};

using rat = basic_rat<int>;