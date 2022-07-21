#include<iostream>
#include<NTL/ZZ.h>
#include"SM3.h"
using namespace NTL;
using namespace std;
ZZ mod, n;
ZZ a, b;     //Elliptic Curve Function: y ^ 2 = x ^ 3 + a * x + b
const int SIZE = 256;
const ZZ pow(uint32_t n) {
	if (n == 0)
		return (ZZ)1;
	return 2 * pow(n - 1);
}
const ZZ modPow(const ZZ& num,ZZ t,const ZZ& n) {
	ZZ temp, ans;
	ans = 1;
	temp = num;
	while (t > 0) {
		if (t % 2 == 1)
			ans = ans * temp % n;
		temp *= temp;
		temp %= n;
		t /= 2;
	}
	return ans;
}
const uint32_t cut(const ZZ& num) {     //Take Low 32 Bits Of ZZ, Transfer Into Uint32_t
	uint32_t n = 0;
	uint32_t temp = 1;
	ZZ x = num;
	for (int i = 0; i < 32; i++) {
		if (x % 2 == 1)
			n += temp;
		temp += temp;
		x /= 2;
	}
	return n;
}
ZZ reverse(const ZZ& num, const ZZ& n) { return modPow(num, n - 2, n); }
class ellPoint {
public:
	ZZ x,y;
	ellPoint(ZZ _x = (ZZ)0, ZZ _y = (ZZ)0) { x = _x; y = _y;}
	const bool setEllPoint(const ZZ& _x) {
		x = _x;
		ZZ num = (x * x * x + a * x + b) % mod;
		if (modPow(num, (mod - 1) / 2, mod) == 1) {
			y = modPow(num, (mod + 1) / 4, mod);
			return 1;
		}
		return 0;
	}
	void pointHash(uint32_t* hash) {     //x,y are 256 bits
		uint32_t XY[16];
		for (uint32_t i = 0; i < 8; i++) {
			XY[i] = cut(x / pow(7 - i));
			XY[i + 8] = cut(y / pow(7 - i));
		}
	}
	const bool OnCurve() { return y * y % mod == (x * x * x + a * x + b) % mod; }
	const ellPoint operator+(const ellPoint& n)const {
		ellPoint p;
		ZZ lambda;
		if (x == n.x && y == n.y)
			lambda = (3 * x * x + a) * reverse(2 * y, mod) % mod;
		else
			lambda = (n.y - y) * reverse(n.x - x, mod) % mod;
		p.x = (lambda * lambda - x - n.x) % mod;
		p.y = (lambda * (x - p.x) - y) % mod;
		return p;
	}
	const ellPoint& operator+=(const ellPoint& n) { *this = *this + n; return *this; }
	const ellPoint operator*(const ZZ& n)const {
		ellPoint p;
		ellPoint temp = *this;
		ZZ num = n;
		while (num % 2 == 0) {
			num /= (ZZ)2;
			temp += temp;
		}
		p = temp;
		while (num > 0) {
			num /= 2;
			temp += temp;
			if (num % 2 == 1)
				p += temp;
		}
		return p;
	}
	friend const ellPoint operator*(const ZZ& n, const ellPoint& p) { return p * n; }
	const ellPoint& operator*=(const ZZ& n) { *this = *this * n; return *this; }
	friend ostream& operator<<(ostream& os, const ellPoint& n) {
		os << "(" << n.x << ", " << n.y << ")";
		return os;
	}
	const bool operator!=(const ellPoint& p) { if (x == p.x && y == p.y) return 0; return 1; }
	void place(uint32_t* pos)const {
		ellPoint p = *this;
		uint32_t* mem = pos;
		for (int i = SIZE / 32 - 1; i >= 0; i--) {
			mem[i] = cut(p.x);
			p.x /= pow((uint32_t)32);
		}
		mem = &pos[SIZE / 32];
		for (int i = SIZE / 32 - 1; i >= 0; i--) {
			mem[i] = cut(p.y);
			p.y /= pow((uint32_t)32);
		}
	}
};
const ZZ transfer(const uint32_t* ptr) {
	ZZ ans = (ZZ)0;
	ZZ temp;
	for (int i = 0; i < 8; i++) {
		if (ptr[i] >= 2147483648) {
			temp = ptr[i];
			temp += pow((uint32_t)32);
		}
		else
			temp = ptr[i];
		ans += temp * pow((uint32_t)(7 - i) * 32);
	}
	return ans;
}
static struct RS {
	ZZ r, s;
};
const RS Sign(const uint32_t* Z_AM_1,const ZZ& d_A,const ZZ& Xa,const ZZ& k,const ZZ& n) {
	uint32_t e[8];
	sm3_context ctx(72);
	ctx.sm3_hash(Z_AM_1, e);
	RS sig;
	sig.r = (transfer(e) + Xa) % n;
	sig.s = (reverse(1 + d_A, n) * (k - sig.r * d_A)) % n;
	return sig;
}
const bool Verify(const RS& sig, const ellPoint& G,const ZZ& k, const uint32_t* Z_AM,const uint64_t& len,const ZZ& n) {
	ZZ deduceK = (k - sig.s) * reverse(sig.s + sig.r, n) % n;
	ellPoint P_A = deduceK * G;
	sm3_context ctx(len);
	uint32_t e[8];
	ctx.sm3_hash(Z_AM, e);
	ellPoint temp = sig.s * G + (sig.r + sig.s) % n * P_A;
	if (sig.r%n == (transfer(e) + temp.x) % n)
		return 1;
	return 0;
}
void sm2_test(const uint32_t* arra,const uint32_t* arrb) {
	const uint32_t arrn[8] = { 0x8542D69E, 0x4C044F18, 0xE8B92435, 0xBF6FF7DD, 0x29772063, 0x0485628D, 0x5AE74EE7, 0xC32E79B7 };
	ZZ n = transfer(arrn);
	ellPoint G, temp;
	const uint32_t arrX_G[8] = { 0x421DEBD6, 0x1B62EAB6, 0x746434EB, 0xC3CC315E, 0x32220B3B, 0xADD50BDC, 0x4C4E6C14, 0x7FEDD43D };
	const uint32_t arrY_G[8] = { 0x0680512B, 0xCBB42C07, 0xD47349D2, 0x153B70C4, 0xE5D7FDFC, 0xBFA36EA1, 0xA85841B9, 0xE46E09A2 };
	G.x = transfer(arrX_G); G.y = transfer(arrY_G);
	cout << "n = " << n << endl;
	cout << "G = " << G << endl;
	uint32_t ID_A = (uint32_t)"ID_A";
	uint32_t ENTL_A = 32;
	sm3_context ctx(32 * 4 + SIZE * 4), _ctx(72);
	uint32_t message[2 + 2 * 8 + SIZE / 32 * 4] = { ENTL_A,ID_A};
	for (int i = 2; i < 10; i++)
		message[i] = arra[i - 2];
	for (int i = 8; i < 18; i++)
		message[i] = arrb[i - 8];
	G.place(&message[18]);
	ZZ k = RandomBits_ZZ(160);
	cout << "Random k Select: " << k << endl;
	ellPoint XYa = k * G;
	XYa.place(&message[18 + SIZE / 32 * 2]);
	uint32_t Z_AM_1[9];
	ctx.sm3_hash(message, Z_AM_1);
	Z_AM_1[8] = (int)" :) ";//A Simple Message
	ZZ d_A = (ZZ)12345;     //Random Key
	RS sig = Sign(Z_AM_1, d_A, XYa.x, k, n);
	if (Verify(sig, G, k, Z_AM_1, 72, n))
		cout << "Check!" << endl;
	else
		cout << "Incorrect!" << endl;
}
//Assume that every message in set is 256 bits
void EllipticCurveMultisetHash(uint32_t** messageSet, const uint32_t& numOfElements, uint32_t* hash) {     //https://arxiv.org/pdf/1601.06502v1.pdf
	for (int i = 0; i < 8; i++)
		hash[i] = 0;
	sm3_context ctx(256),_ctx(288);
	uint32_t output[9],_output[8];
	ellPoint temp;
	for (int i = 0; i < 8; i++)
		hash[i] = 0;
	for (uint32_t i = 0; i < numOfElements; i++) {
		ctx.sm3_hash(messageSet[i], &output[1]);
		for (output[0] = 0; output[0] < pow(31); output[0]++) {
			_ctx.sm3_hash(output, _output);
			if(temp.setEllPoint(transfer(_output)))     //check if this x is on the elliptic curve
				break;
		}
		temp.pointHash(_output);
		for (int i = 0; i < 8; i++)
			hash[i] ^= _output[i];
	}
}
inline bool check(const uint32_t* h1, const uint32_t* h2) {
	for (int i = 0; i < 8; i++)
		if (h1[i] != h2[i])
			return 0;
	return 1;
}
#include<cstdlib>
void ECMHtest() {     //Test The Homomorphic Of ECMH Function
	uint32_t** Message = new uint32_t*[3];
	for (int i = 0; i < 3; i++)
		Message[i] = new uint32_t[8];
	uint32_t hash123[8];
	uint32_t hash1[8], hash12[8], hash23[8], hash3[8];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			Message[i][j] = rand();     //filled with random message
	EllipticCurveMultisetHash(Message, 3, hash123);
	EllipticCurveMultisetHash(Message, 1, hash1);
	EllipticCurveMultisetHash(Message, 2, hash12);
	EllipticCurveMultisetHash(&Message[1], 2, hash23);
	EllipticCurveMultisetHash(&Message[2], 1, hash3);
	for (int i = 0; i < 3; i++)
		delete[] Message[i];
	delete[] Message;
	for (int i = 0; i < 8; i++) {     //compute hash
		hash1[i] ^= hash23[i];
		hash12[i] ^= hash3[i];
	}
	if (check(hash1, hash123) && check(hash12, hash123)) {
		cout << "Check!" << endl;
		return;
	}
	cout << "Incorrect!";
}
int main()
{
	const uint32_t arrmod[8] = { 0x8542D69E, 0x4C044F18 ,0xE8B92435 ,0xBF6FF7DE ,0x45728391 ,0x5C45517D ,0x722EDB8B ,0x08F1DFC3 };
	mod = transfer(arrmod);
	const uint32_t arra[8] = { 0x787968B4 ,0xFA32C3FD ,0x2417842E ,0x73BBFEFF ,0x2F3C848B ,0x6831D7E0 ,0xEC65228B ,0x3937E498 };
	const uint32_t arrb[8] = { 0x63E4C6D3 ,0xB23B0C84 ,0x9CF84241 ,0x484BFE48 ,0xF61D59A5 ,0xB16BA06E ,0x6E12D1DA ,0x27C5249A };
	a = transfer(arra); b = transfer(arrb);
	//sm2_test(arra, arrb);
	//ECMHtest();
	return 0;
}
