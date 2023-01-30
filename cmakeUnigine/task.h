// THIS FILE IS STUB FOR DECLARATION FUCTION
// DO NOT MODIFICATE IT - IT`S PART OF ANOTHER PROJECT
#pragma once
#include "common.h"
#include <vector>
#include <iostream>
namespace Task
{
void checkVisible(const std::vector<unit> &input_units, std::vector<int> &result);
}
//#include "task.h"
//#include <cmath>
//#include <thread>
//#include <xmmintrin.h>
//#include <pmmintrin.h>
//#include <immintrin.h>
//const float PI = 3.14159f;
//const float rad = 3.14159f / 180.0f;
//const float half = 1.0f / 2.0f;
//void rotateSSE(__m128 &borderVecs,const vec2& vec, const float& angle)
//{
//	float cosine = cosf(angle * rad);
//	float sinus = sinf(angle * rad);
//	__m128 border = _mm_setr_ps(vec.x, vec.y, vec.x, vec.y);
//	__m128 rotateMx = _mm_setr_ps(cosine, -sinus, sinus, cosine);
//	__m128 rotated1 = _mm_mul_ps(border, rotateMx);
//	rotateMx = _mm_shuffle_ps(rotateMx, rotateMx, 0x18);//latency 1
//	__m128 rotated2 = _mm_mul_ps(border,rotateMx);
//	//borderVecs = _mm_hadd_ps(rotated1, rotated2);//поменять hadd 7 latency CPI 2
//	borderVecs = _mm_add_ps(_mm_permute_ps(rotated1,_MM_SHUFFLE(0,3,1,2)), _mm_permute_ps(rotated2,_MM_SHUFFLE(1,2,0,3)));//add latency 4 CPI 0.5 permute latency 1 x2 CPI 1
//	//borderVecs = _mm_shuffle_ps(borderVecs, borderVecs, 0xB1);
//}
//bool checkInCircleSSE(__m128 sub, float& squareRadius)
//{
//	__m128 mulVecs = _mm_mul_ps(sub, sub);
//	__m128 addVecs = _mm_add_ps(mulVecs, _mm_permute_ps(mulVecs,0xB1));
//	return _mm_comile_ss(addVecs, _mm_set_ss(squareRadius));
//}
//bool checkSectorSSE(__m128 &border, __m128 &vec)
//{
//	__m128 mulVecs = _mm_mul_ps(border, vec);
//	__m128 subVecs = _mm_sub_ps(mulVecs, _mm_movehdup_ps(mulVecs));
//	__m128 cmp = _mm_cmplt_ps(subVecs, _mm_setzero_ps());
//	return _mm_movemask_ps(cmp) == 4;
//}
//
//__m128 substractSSE(const vec2& v1, const vec2& v2)
//{
//	return _mm_sub_ps(_mm_castpd_ps(_mm_load1_pd((double*)&v1)), _mm_castpd_ps(_mm_load1_pd((double*)&v2)));
//}
//auto checkUnitVision(const std::vector<unit>& input_units, std::vector<int>& result, int startFrom, int step)
//{
//	int size = input_units.size();
//	for (int i = startFrom; i < size; i += step)
//	{
//		__m128 borderVecs128;
//		rotateSSE(borderVecs128, input_units[i].direction, input_units[i].fov_deg * half);
//		float squareRadius = input_units[i].distance * input_units[i].distance;
//		for (int j = 0; j < size; ++j)
//		{
//			__m128 sub = substractSSE(input_units[j].position, input_units[i].position);
//			if (checkInCircleSSE(sub, squareRadius))
//			{
//				if (checkSectorSSE(borderVecs128, sub))
//				{
//					//std::cout << i <<"sees" << j<<'\n';
//					++result[i];
//				}
//			}
//		}
//	}
//}
//void Task::checkVisible(const std::vector<unit>& input_units, std::vector<int>& result)
//{
//	result.resize(input_units.size());
//	int thread_num = std::thread::hardware_concurrency();
//	std::vector<std::thread> threads(thread_num-1);
//	for (int threadIdx = 0; threadIdx < thread_num - 1; ++threadIdx)
//		threads[threadIdx] = std::thread(checkUnitVision, std::ref(input_units), std::ref(result), threadIdx, thread_num);
//	checkUnitVision(input_units, result, thread_num - 1, thread_num);
//	//checkUnitVision(input_units, result, 0, 1); 
//	for (auto& t : threads)
//		if(t.joinable())
//		t.join();
//}