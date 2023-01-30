// THIS FILE IS ONLY FOR STRUCT DECLARATION
// DO NOT MODIFICATE IT - IT`S PART OF ANOTHER PROJECT

// это файл только для деклараций общих структур
// не модифицируйте его - возможно это часть другого проекта

#pragma once

struct vec2
{
	float x = 0.0f;
	float y = 0.0f;
};

struct unit
{
	vec2 position; // position of unit (-10^5...10^5, -10^5...10^5)
	vec2 direction; // normalized view direction
	float fov_deg = 0.0f; // horizontal field-of-view in degrees (0...180)
	float distance = 0.0f; // view distance (0...10^5)
};
//#include "task.h"
//#include <cmath>
//#include <thread>
//#include <iterator>
//
//#include <xmmintrin.h>
//#include <pmmintrin.h>
//#include <immintrin.h>
//
//const float radd_coeff = 0.0174533 / 2;
//
//void rotate(__m256& border_vectors, __m256& direction_vector, const float& turn_angle1, const float& turn_angle2)
//{
//	//float sincos[8] = { std::cos(turn_angle1),std::sin(turn_angle1), };
//	float cosA1 = std::cos(turn_angle1);
//	float sinA1 = std::sin(turn_angle1);
//	float cosA2 = std::cos(turn_angle2);
//	float sinA2 = std::sin(turn_angle2);
//
//	__m256 rotate_vec = _mm256_setr_ps(cosA1, -sinA1, sinA1, cosA1, cosA2, -sinA2, sinA2, cosA2);
//
//	__m256 rotate1 = _mm256_mul_ps(direction_vector, rotate_vec);
//	rotate_vec = _mm256_permute_ps(rotate_vec, _MM_SHUFFLE(3, 1, 2, 0));
//	__m256 rotate2 = _mm256_mul_ps(direction_vector, rotate_vec);
//
//	border_vectors = _mm256_add_ps(_mm256_permute_ps(rotate1, _MM_SHUFFLE(0, 3, 1, 2)), _mm256_permute_ps(rotate2, _MM_SHUFFLE(1, 2, 0, 3)));
//}
//
//int unit_in_circle(__m256& unit_vec, __m256& R2)
//{
//	__m256 distance2 = _mm256_mul_ps(unit_vec, unit_vec);
//	distance2 = _mm256_add_ps(distance2, _mm256_permute_ps(distance2, _MM_SHUFFLE(2, 3, 0, 1)));
//	__m256 cmp = _mm256_cmp_ps(distance2, R2, 2);
//	return _mm256_movemask_ps(cmp);
//}
//
//int unit_in_sector(__m256& border_vectors, __m256& unit_vec)
//{
//	__m256 mul = _mm256_mul_ps(border_vectors, unit_vec);
//	__m256 shuf = _mm256_movehdup_ps(mul);
//	__m256 sub = _mm256_sub_ps(mul, shuf);
//	__m256 cmp = _mm256_cmp_ps(sub, _mm256_setzero_ps(), 1);
//	return _mm256_movemask_ps(cmp);
//}
//
//void check_units(const std::vector<unit>& input_units, std::vector<int>& result, int thrd_ix, int step)
//{
//	int units_count = input_units.size();
//	//int last_obs = 0;
//	//if (units_count % 2)
//	// last_obs = —units_count;
//
//	for (int obs_ix = thrd_ix * 2; obs_ix < units_count; obs_ix += step)
//	{
//		__m256 observer_position = _mm256_setr_m128(_mm_castpd_ps(_mm_load1_pd((double*)&input_units[obs_ix].position)), _mm_castpd_ps(_mm_load1_pd((double*)&input_units[obs_ix + 1].position)));
//		__m256 direction_vector = _mm256_setr_m128(_mm_castpd_ps(_mm_load1_pd((double*)&input_units[obs_ix].direction)), _mm_castpd_ps(_mm_load1_pd((double*)&input_units[obs_ix + 1].direction)));
//
//		__m256 border_vectors;
//		rotate(border_vectors, direction_vector, input_units[obs_ix].fov_deg * radd_coeff, input_units[obs_ix + 1].fov_deg * radd_coeff);
//
//		__m256 R2 = _mm256_setr_m128(_mm_load_ss((float*)&input_units[obs_ix].distance), _mm_load_ss((float*)&input_units[obs_ix + 1].distance));
//		R2 = _mm256_mul_ps(R2, R2);
//
//		for (auto& unit : input_units)
//		{
//			__m128 unit_position = _mm_castpd_ps(_mm_load1_pd((double*)&unit.position));
//			__m256 unit_vector = _mm256_set_m128(unit_position, unit_position);
//			unit_vector = _mm256_sub_ps(unit_vector, observer_position);
//			int in_circle_mask = unit_in_circle(unit_vector, R2);
//			if (in_circle_mask)
//			{
//				int in_sector_mask = unit_in_sector(border_vectors, unit_vector);
//				if ((in_circle_mask & 1) && (in_sector_mask & 15) == 4)
//					++result[obs_ix];
//				if ((in_circle_mask & 16) && (in_sector_mask & 240) == 64)
//					++result[obs_ix + 1];
//			}
//		}
//	}
//
//}
//
//void Task::checkVisible(const std::vector<unit>& input_units, std::vector<int>& result)
//{
//	int units_count = input_units.size();
//	result.resize(units_count, 0);
//
//	int threads_count = 3;
//
//	int step = (threads_count + 1) * 2;
//
//	std::vector<std::thread> threads(threads_count);
//
//	int thrd_ix = 0;
//	for (auto& thread : threads)
//	{
//		thread = std::thread(check_units, std::ref(input_units), std::ref(result), thrd_ix, step);
//		++thrd_ix;
//	}
//
//	check_units(input_units, result, thrd_ix, step);
//
//	for (auto& thread : threads)
//		thread.join();
//}