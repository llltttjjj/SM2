/*This file is used as hash function in sm2.cpp*/
#include<iostream>
using namespace std;
static const uint64_t pow(int n) {
	if (n == 0)
		return 1;
	return 2 * pow(n - 1);
}
class sm3_context {
private:
	uint32_t* info;     //unknown bits message
	uint64_t size;     //size of message
	unsigned int messlen;
	uint32_t dig[8];     //256 bits digest
	const uint64_t U32 = pow(32);
	void sm3_init(const uint32_t* message) {
		uint64_t tmp = (size - 1) / 32 + 1;
		for (uint64_t i = 0; i < tmp; i++)
			info[i] = message[i];
		for (uint64_t i = tmp; i < messlen; i++)
			info[i] = 0;
		if (size % 32 != 0)     //在消息最后一位置一
			info[tmp - 1] += pow(int(31 - size % 32));
		else
			info[tmp] += U32 / 2;
		info[messlen - 2] = size / U32;
		info[messlen - 1] = size % U32;
	}
	uint32_t ls(const uint32_t val, int num) { return val >> (32 - num) | (val << num); }
	uint32_t P0(const uint32_t val) { return val ^ ls(val, 9) ^ ls(val, 17); }
	uint32_t P1(const uint32_t val) { return val ^ ls(val, 15) ^ ls(val, 23); }
	uint32_t T(short j) { if (j < 16) return 0x79CC4519; return 0x7A879D8A; }
	uint32_t FF(uint32_t X, uint32_t Y, uint32_t Z, short j) {
		if (j < 16) return X ^ Y ^ Z;
		return (X & Y) | (X & Z) | (Y & Z);
	}
	uint32_t GG(uint32_t X, uint32_t Y, uint32_t Z, short j) {
		if (j < 16) return X ^ Y ^ Z;
		return (X & Y) | (~X & Z);
	}
	//message is 512 bits, IV is 256 bits
	void itercompress(const uint32_t* message, const uint32_t* IV, uint32_t* output) {
		uint32_t W[68];
		uint32_t W_[64];
		for (int j = 0; j < 16; j++)
			W[j] = message[j];
		for (int j = 16; j < 68; j++)
			W[j] = P1(W[j - 16] ^ W[j - 9] ^ ls(W[j - 3], 15)) ^ ls(W[j - 13], 7) ^ W[j - 6];
		for (int j = 0; j < 64; j++)
			W_[j] = W[j] ^ W[j + 4];
		uint32_t SS1, SS2, TT1, TT2;
		for (int i = 0; i < 8; i++)
			output[i] = IV[i];
		uint32_t& A = output[0];
		uint32_t& B = output[1];
		uint32_t& C = output[2];
		uint32_t& D = output[3];
		uint32_t& E = output[4];
		uint32_t& F = output[5];
		uint32_t& G = output[6];
		uint32_t& H = output[7];
		for (int j = 0; j < 64; j++) {
			SS1 = ls(ls(A, 12) + E + ls(T(j), j), 7);
			SS2 = SS1 ^ ls(A, 12);
			TT1 = FF(A, B, C, j) + D + SS2 + W_[j];
			TT2 = GG(E, F, G, j) + H + SS1 + W[j];
			D = C;
			C = ls(B, 9);
			B = A;
			A = TT1;
			H = G;
			G = ls(F, 19);
			F = E;
			E = P0(TT2);
		}
		for (int i = 0; i < 8; i++)
			output[i] = IV[i] ^ output[i];
	}
	void reset() {
		dig[0] = 0x7380166F;
		dig[1] = 0x4914B2B9;
		dig[2] = 0x172442D7;
		dig[3] = 0xDA8A0600;

		dig[4] = 0xA96F30BC;
		dig[5] = 0x163138AA;
		dig[6] = 0xE38DEE4D;
		dig[7] = 0xB0FB0E4E;
	}
public:
	sm3_context(const uint64_t sizeOfMessage) {
		messlen = ((sizeOfMessage - 447) / 512 + 2) * 16;
		info = new uint32_t[messlen];
		size = sizeOfMessage;
		reset();
	}
	void sm3_extendcontext(const uint64_t sizeOfMessage, const uint32_t* hash) {
		messlen = ((sizeOfMessage - 447) / 512 + 2) * 16;
		info = new uint32_t[messlen];
		size = sizeOfMessage;
		for (int i = 0; i < 8; i++)
			dig[i] = hash[i];
	}
	void sm3_hash(const uint32_t* message, uint32_t* output) {
		sm3_init(message);
		for (int i = 0; i < messlen / 16; i++) {
			itercompress(&info[i * 16], dig, output);
			for (int j = 0; j < 8; j++)
				dig[i] = output[i];
		}
		reset();
	}
	~sm3_context() {
		delete[] info;
	}
};
