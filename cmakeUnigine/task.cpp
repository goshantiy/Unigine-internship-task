#include "task.h"
#include <cmath>
#include <thread>
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <immintrin.h>
const float PI = 3.14159f;
const float rad = 3.14159f / 180.0f;
const float half = 1.0f / 2.0f;
static const unsigned int grid_size = 16;
__m128 grid_sizeSSE = _mm_set1_ps(grid_size);
const float coeff = 1.3f;

struct Cell
{
	std::vector<int> units;
};
struct Grid
{
	Cell grid[grid_size+1][grid_size+1];
	float stepX;
	float stepY;
	float xBorder;
	float yBorder;
	__m128 gridBounds;
	__m128 gridStep;
};
void rotateSSE(__m128& borderVecs, const vec2& vec, const float& angle)
{
	float cosine = cosf(angle * rad);
	float sinus = sinf(angle * rad);
	__m128 border = _mm_castpd_ps(_mm_load1_pd((double*)&vec));
	__m128 rotateMx = _mm_setr_ps(cosine, -sinus, sinus, cosine);
	__m128 rotated1 = _mm_mul_ps(border, rotateMx);
	rotateMx = _mm_permute_ps(rotateMx, _MM_SHUFFLE(3, 1, 2, 0));
	__m128 rotated2 = _mm_mul_ps(border, rotateMx);
	borderVecs = _mm_add_ps(_mm_permute_ps(rotated1, _MM_SHUFFLE(0, 3, 1, 2)), _mm_permute_ps(rotated2, _MM_SHUFFLE(1, 2, 0, 3)));
}
bool checkInCircleSSE(const __m128& sub, const __m128& squareRadius)
{
	__m128 mulVecs = _mm_mul_ps(sub, sub);
	__m128 addVecs = _mm_add_ps(mulVecs, _mm_permute_ps(mulVecs,_MM_SHUFFLE(2,3,0,1)));
	return _mm_comile_ss(addVecs, squareRadius);
}
bool checkSectorSSE(const __m128& border, const __m128& vec)
{
	__m128 mulVecs = _mm_mul_ps(border, vec);
	__m128 subVecs = _mm_sub_ps(mulVecs, _mm_movehdup_ps(mulVecs));
	__m128 cmp = _mm_cmplt_ps(subVecs, _mm_setzero_ps());
	return _mm_movemask_ps(cmp) == 4;
}
__m128 substractSSE(const __m128& v1, const __m128& v2)
{
	return _mm_sub_ps(v1,v2);
}
void fillGrid(Grid& grid, const std::vector<unit>& input_units)
{
	float xmin = std::numeric_limits<float>::max();
	float xmax = std::numeric_limits<float>::lowest();
	float ymin = std::numeric_limits<float>::max();
	float ymax = std::numeric_limits<float>::lowest();

	for (int i = 0; i < input_units.size(); ++i)
	{
		if (input_units[i].position.x < xmin)
			xmin = input_units[i].position.x;
		if (input_units[i].position.x > xmax)
			xmax = input_units[i].position.x;
		if (input_units[i].position.y < ymin)
			ymin = input_units[i].position.y;
		if (input_units[i].position.y > ymax)
			ymax = input_units[i].position.y;
	}
	grid.stepX = (xmax - xmin) / (grid_size);
	grid.stepY = (ymax - ymin) / (grid_size);
	grid.xBorder = xmin;
	grid.yBorder = ymin;
	for (int i = 0; i < input_units.size(); ++i)
	{
		unsigned int xpos = ((input_units[i].position.x) - (xmin)) / grid.stepX;
		unsigned int ypos = ((input_units[i].position.y) - (ymin)) / grid.stepY;
		grid.grid[xpos][ypos].units.push_back(i);
	}
}

auto checkUnitVision(const std::vector<unit>& input_units, std::vector<int>& result, int startFrom, int step, std::vector<__m128>& xmmpos, Grid& grd)
{
	int size = input_units.size();
	for (int i = startFrom; i < size; i += step)
	{
		__m128 borderVecs128;
		rotateSSE(borderVecs128, input_units[i].direction, input_units[i].fov_deg / 2);
		//__m128 dist = _mm_setr_ps(input_units[i].distance, input_units[i].distance,-input_units[i].distance, -input_units[i].distance);
		//__m128 bbPoints = _mm_sub_ps(xmmpos[i], dist);//yxyx - dst
		//__m128 sub = _mm_sub_ps(bbPoints, grd.gridBounds);//yxyx - gridBords
		//__m128 div = _mm_div_ps(sub, grd.gridStep);//yxyx / step
		////__m128 cmplt = _mm_cmplt_ps(div, _mm_setzero_ps());//yxyx < 0 ? 1 : 0
		////__m128 blendlt = _mm_blendv_ps(div, _mm_setzero_ps(), cmplt);//if 1 than 0 else no change
		////__m128 cmpgt = _mm_cmpgt_ps(blendlt, grid_sizeSSE);//yxyx > gridsize ? 1:0
		////__m128 blendgt = _mm_blendv_ps( blendlt, grid_sizeSSE, cmpgt);//if 1 gridSize else no change
		//__m128 max = _mm_max_ps(div, _mm_setzero_ps());
		//__m128 min = _mm_min_ps(max, grid_sizeSSE);
		//float yxyx[4];
		//_mm_store_ps(yxyx,min);
		float xmin, xmax, ymin, ymax;

		float borders[4];
		_mm_store_ps(borders, borderVecs128);
		xmin = input_units[i].position.x;
		xmax = input_units[i].position.x;
		ymin = input_units[i].position.y;
		ymax = input_units[i].position.y;

		vec2 vr;
		vr.x = (borders[1] * input_units[i].distance * coeff) + input_units[i].position.x;
		vr.y = (borders[0] * input_units[i].distance * coeff) + input_units[i].position.y;
			
		vec2 vl;
		vl.x = (borders[3] * input_units[i].distance * coeff) + input_units[i].position.x;
		vl.y = (borders[2] * input_units[i].distance * coeff) + input_units[i].position.y;

		vec2 vd = input_units[i].direction;
		vd.x = (vd.x * input_units[i].distance * coeff) + input_units[i].position.x;
		vd.y = (vd.y * input_units[i].distance * coeff) + input_units[i].position.y;

		if (vr.x > xmax)
			xmax = vr.x;
		if (vr.x < xmin)
			xmin = vr.x;
		if (vr.y > ymax)
			ymax = vr.y;
		if (vr.y < ymin)
			ymin = vr.y;

		if (vl.x > xmax)
			xmax = vl.x;
		if (vl.x < xmin)
			xmin = vl.x;
		if (vl.y > ymax)
			ymax = vl.y;
		if (vl.y < ymin)
			ymin = vl.y;

		if (vd.x > xmax)
			xmax = vd.x;
		if (vd.x < xmin)
			xmin = vd.x;
		if (vd.y > ymax)
			ymax = vd.y;
		if (vd.y < ymin)
			ymin = vd.y;
		int xposmin = (xmin - grd.xBorder) / grd.stepX;
		if (xposmin < 0)
			xposmin = 0;
		int xposmax = (xmax - grd.xBorder) / grd.stepX;
		if (xposmax > grid_size)
			xposmax = grid_size;

		int yposmin = (ymin - grd.yBorder) / grd.stepY;
		if (yposmin < 0)
			yposmin = 0;
		int yposmax = (ymax - grd.yBorder) / grd.stepY;
		if (yposmax > grid_size)
			yposmax = grid_size;
		float squareRadius = input_units[i].distance * input_units[i].distance;
		__m128 radius2 = _mm_set1_ps(squareRadius);
		for (int ix = xposmin; ix <= xposmax; ++ix)
			for (int jx = yposmin; jx <= yposmax; ++jx)
			{
				for (int sz = 0; sz < grd.grid[ix][jx].units.size(); ++sz) 
				{
					__m128 sub = substractSSE(xmmpos[grd.grid[ix][jx].units[sz]], xmmpos[i]);
						if (checkInCircleSSE(sub, radius2))
						{
							if (checkSectorSSE(borderVecs128,sub))
							{
								++result[i];
							}
						}
				}
			}

	}
}
void castVec2ToXmm(std::vector<__m128>& dst, vec2& src)
{
	dst.push_back(_mm_castpd_ps(_mm_load1_pd((double*)&src)));
}
void Task::checkVisible(const std::vector<unit>& input_units, std::vector<int>& result)
{
	Grid grid;
	fillGrid(grid, input_units);//threads?
	grid.gridBounds = _mm_castpd_ps(_mm_load_pd1((double*)&grid.xBorder));//struct linear memory
	grid.gridStep = _mm_castpd_ps(_mm_load_pd1((double*)&grid.stepX));//struct linear memory
	result.resize(input_units.size());
	int thread_num = std::thread::hardware_concurrency();
	std::vector<__m128> xmmpos;
	for (auto it : input_units)
		castVec2ToXmm(xmmpos, it.position);
	std::vector<std::thread> threads(thread_num - 1);
	for (int threadIdx = 0; threadIdx < thread_num - 1; ++threadIdx)
		threads[threadIdx] = std::thread(checkUnitVision, std::ref(input_units), std::ref(result), threadIdx, thread_num, std::ref(xmmpos), std::ref(grid));
	checkUnitVision(input_units, result, thread_num - 1, thread_num, xmmpos, grid);
	//checkUnitVision(input_units, result, 0, 1, xmmpos, grid);
	for (auto& t : threads)
		if (t.joinable())
			t.join();
}