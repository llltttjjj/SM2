#ifndef GOOGLEPASSWORDCHECKUP_H_
#define GOOGLEPASSWORDCHECKUP_H_



#include<iostream>
#include<NTL/ZZ.h>
#include<cstdlib>
#include"SM3.h"
using namespace NTL;
using namespace std;
static const ZZ pow(uint32_t n) {
	if (n == 0)
		return (ZZ)1;
	return 2 * pow(n - 1);
}
static const ZZ modPow(const ZZ& num, ZZ t, const ZZ& n) {
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
/* This checkup takes 2 Bytes of h_i, 1M random username and key, 2^4 pairs of h^b in set S. */
static const int Size = 1048576;     // 1M
uint32_t randomU_K[Size][2];     //[0] for username, [1] for key
struct k_v{
	uint16_t hashRec;
	ZZ v;
};
k_v memory[Size];
const ZZ Pow(const ZZ& h, ZZ b) {
	ZZ temp, ans;
	ans = 1;
	temp = h;
	while (b > 0) {
		if (b % 2 == 1)
			ans = ans * temp;
		temp *= temp;
		b /= 2;
	}
	return ans;
}
static const ZZ ProcessDataInfo() {
	cout << "Random Fill: " << endl;
	for (int i = 0; i < Size; i++) {
		randomU_K[i][0] = rand();
		randomU_K[i][1] = rand();
		cout << (float)i / Size * 100 << "\r";
	}
	cout << endl;
	ZZ b = RandomBits_ZZ(8);
	sm3_context ctx(64);
	cout << "Hash Process: " << endl;
	for (int i = 0; i < Size; i++) {
		static uint32_t output[8];
		ctx.sm3_hash(randomU_K[i], output);
		memory[i].hashRec = output[0] / pow(16);
		memory[i].v = Pow(transfer(output), b);
		cout << (float)i / Size * 100 << "\r";
	}
	cout << endl;
	return b;
}

const ZZ commonFactor(const ZZ& m, const ZZ& n) {
	if (m % n == 0)
		return n;
	if (n == 1)
		return (ZZ)1;
	return commonFactor(n, m % n);
}
static const ZZ reverse(const ZZ& num, const ZZ& n) { return modPow(num, n - 2, n); }
const k_v UserInputNameAndPassword(const uint32_t* U_K,ZZ& N,ZZ& d) {
	ZZ p = RandomPrime_ZZ(1024);
	ZZ q = RandomPrime_ZZ(1024);
	N = p * q;
	ZZ phi = (p - 1) * (q - 1);
	ZZ e;
	do {
		e = 2 * RandomBits_ZZ(256) + 1;
	} while (commonFactor(phi,e)!=1);
	d = InvMod(e, phi);
	sm3_context ctx(64);
	uint32_t output[8];
	ctx.sm3_hash(U_K, output);
	k_v myK_V;
	myK_V.hashRec = output[0] / pow(16);
	myK_V.v = modPow(transfer(output), e, N);
	return myK_V;
}
//return h_ab
const ZZ  FindTheDataSet(const k_v& kv, const ZZ& b, ZZ* VSet,int& num) {
	ZZ h_ab = Pow(kv.v, b);
	num = 0;
	for (int i = 0; i < Size; i++)
		if (memory[i].hashRec == kv.hashRec)
			num++;
	VSet = new ZZ[num];
	int k = 0;
	for (int i = 0; i < Size; i++)
		if (memory[i].hashRec == kv.hashRec) {
			VSet[k] = memory[i].v;
		}
	return h_ab;
}
const bool UserNameAndPasswordDetection(const ZZ& h_ab, const ZZ& N, const ZZ& d, const ZZ* VSet, const int& num) {
	ZZ h_b = modPow(h_ab % N, d, N);
	for (int i = 0; i < num; i++)
		if (h_b == VSet[i])
			return 1;
	return 0;
}
void checkUpTest() {
	ZZ b = ProcessDataInfo();
	for (int i = 0; i < 5; i++) {     //5 Users that their passwords are leaked
		ZZ N, d;
		ZZ* VSet = nullptr;
		int num;
		k_v kv = UserInputNameAndPassword(randomU_K[rand() % Size], N, d);
		ZZ h_ab = FindTheDataSet(kv, b, VSet, num);
		if (UserNameAndPasswordDetection(h_ab, N, d, VSet, num))
			cout << "Leaking detected!" << endl;
		else
			cout << "Didn't found the leaking." << endl;
	}
	for (int i = 0; i < 5; i++) {     //5 Users that their passwords are not leaked
		ZZ N, d;
		ZZ* VSet = nullptr;
		int num;
		uint32_t U_K[2];
		U_K[0] = Pow((ZZ)rand(), (ZZ)3) % pow(31) + rand();
		U_K[1] = Pow((ZZ)rand(), (ZZ)3) % pow(31) + rand();     //take random username and key, unlikely in the random set.
		k_v kv = UserInputNameAndPassword(U_K, N, d);
		ZZ h_ab = FindTheDataSet(kv, b, VSet, num);
		if (UserNameAndPasswordDetection(h_ab, N, d, VSet, num))
			cout << "False leaking detected!" << endl;
		else
			cout << "Didn't found the leaking." << endl;
	}
}



#endif
