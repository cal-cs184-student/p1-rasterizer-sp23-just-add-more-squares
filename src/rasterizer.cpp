#include "rasterizer.h"

using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)


    sample_buffer[y * width + x] = c;
  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    fill_pixel(sx, sy, color);
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  bool inside_triangle_side(float Xi, float Yi,
                            float dXi, float dYi,
                            float x, float y) {
    return -(x - Xi) * dYi + (y - Yi) * dXi >= 0;
  }

  bool inside_triangle(float X0, float Y0,
                       float X1, float Y1,
                       float X2, float Y2,
                       float x, float y) {
    float dX0 = X1 - X0, dY0 = Y1 - Y0;
    float dX1 = X2 - X1, dY1 = Y2 - Y1;
    float dX2 = X0 - X2, dY2 = Y0 - Y2;
    bool L0 = inside_triangle_side(X0, Y0, dX0, dY0, x, y);
    return L0 == inside_triangle_side(X1, Y1, dX1, dY1, x, y)
        && L0 == inside_triangle_side(X2, Y2, dX2, dY2, x, y);
  }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
//    float min_x = floor(min(x0, min(x1, x2))) + 0.5;
//    float max_x = ceil(max(x0, max(x1, x2))) - 0.5;
//    float min_y = floor(min(y0, min(y1, y2))) + 0.5;
//    float max_y = ceil(max(y0, max(y1, y2))) - 0.5;
//
//    for (float x = min_x; x <= max_x; x++) {
//      for (float y = min_y; y <= max_y; y++) {
//        if (inside_triangle(x0, y0, x1, y1, x2, y2, x, y)) {
//          rasterize_point(x, y, color);
//        }
//      }
//    }

    // TODO: Task 2: Update to implement super-sampled rasterization
    double step = 1.0 / (sqrt(this->sample_rate));
    int min_x = floor(min(x0, min(x1, x2))) - 1;
    int max_x = ceil(max(x0, max(x1, x2))) + 1;
    int min_y = floor(min(y0, min(y1, y2))) - 1;
    int max_y = ceil(max(y0, max(y1, y2))) + 1;

    int iters = round(sqrt(this->sample_rate));
    for (int i = min_x; i <= max_x; i++) {
      for (int j = min_y; j <= max_y; j++) {
        Color c = Color(0, 0, 0);
        bool inside = false;
        for (int x = 0; x < iters; x++) {
          for (int y = 0; y < iters; y++) {
            int u = i + (x * step) + (step / 2);
            int v = j + (y * step) + (step / 2);
            if (inside_triangle(x0, y0, x1, y1, x2, y2, u, v)) {
              inside = true;
              c += color * (1.0 / this->sample_rate);
            } else {
              c += Color::White * (1.0 / this->sample_rate);
            }
          }
        }
        double center_x = i, center_y = j;
        if (inside) {
          rasterize_point(center_x, center_y, c);
        }
      }
    }
  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle



  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle




  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support
    this->sample_rate = rate;
    this->sample_buffer.resize(width * height, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;


    this->sample_buffer.resize(width * height, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support


    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        Color col = sample_buffer[y * width + x];

        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        }
      }
    }

  }

  Rasterizer::~Rasterizer() { }


}// CGL
