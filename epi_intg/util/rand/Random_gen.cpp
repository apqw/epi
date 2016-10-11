/*
 * Random_gen.cpp
 *
 *  Created on: 2016/10/11
 *      Author: yasu7890v
 */

#include "Random_gen.h"
#include <random>
#include <algorithm>
#include <array>
#include <functional>

/**
 *  @file 乱数に関する実装
 */

/**
 *  乱数を生成する
 *  @attention スレッドセーフである。(内部でstatic thread_localな乱数生成器が作成される)
 */
real genrand_real()
{
    /**
    *  @brief 乱数生成クラス
    *
    *  コンストラクタでシード等の初期化をして、
    *  外部公開用関数内にstatic thread_localな乱数発生器のインスタンスを置く目的で作成。
    */
    class Random_gen {
        std::mt19937_64 mt;
        std::uniform_real_distribution<real> dist;

        /**
        *  初期化
        *  @todo:固定シードを与える
        */
        inline void init() {
            std::random_device rng;
            std::array<uint32_t, 16> rng_seed;
            std::generate(rng_seed.begin(), rng_seed.end(), std::ref(rng));
            std::seed_seq rng_seed_seq(rng_seed.begin(), rng_seed.end());
            mt.seed(rng_seed_seq);
        }
    public:
        /**
        *  ctor
        */
        Random_gen() :dist(0.0, 1.0) {
            init();
        }

        /**
        *  乱数を生成する
        */
        double gen_rand_real() {
            return dist(mt);
        }
    };
	static thread_local Random_gen rng;
	return rng.gen_rand_real();
}
