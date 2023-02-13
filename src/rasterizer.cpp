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

    int super_x = x * sqrt(sample_rate);
    int super_y = y * sqrt(sample_rate);
    for (int i = super_x; i < super_x + sqrt(sample_rate); i++) {
      for (int j = super_y; j < super_y + sqrt(sample_rate); j++) {
        sample_buffer[j * width * sqrt(sample_rate) + i] = c;
      }
    }
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
    float step = 1.0 / sqrt(sample_rate);
    float min_x = floor(min(x0, min(x1, x2))) + 0.5 * step;
    float max_x = ceil(max(x0, max(x1, x2))) - 0.5 * step;
    float min_y = floor(min(y0, min(y1, y2))) + 0.5 * step;
    float max_y = ceil(max(y0, max(y1, y2))) - 0.5 * step;

    for (float x = min_x; x <= max_x; x += step) {
      for (float y = min_y; y <= max_y; y += step) {
        if (inside_triangle(x0, y0, x1, y1, x2, y2, x, y)) {
          int i = round((x - 0.5 * step) / step);
          int j = round((y - 0.5 * step) / step);
          sample_buffer[j * width * sqrt(sample_rate) + i] = color;
        }
      }
    }


  }

    Vector3D barycentric(float x0, float y0,
                         float x1, float y1,
                         float x2, float y2,
                         float x, float y) {
      double alpha = (-(x - x1) * (y2 - y1) + (y - y1) * (x2 - x1)) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
      double beta = (-(x - x2) * (y0 - y2) + (y - y2) * (x0 - x2)) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
      double gamma = 1 - alpha - beta;
      return Vector3D(alpha, beta, gamma);
    }

  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2) {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle
    float step = 1.0 / sqrt(sample_rate);
    float min_x = floor(min(x0, min(x1, x2))) + 0.5 * step;
    float max_x = ceil(max(x0, max(x1, x2))) - 0.5 * step;
    float min_y = floor(min(y0, min(y1, y2))) + 0.5 * step;
    float max_y = ceil(max(y0, max(y1, y2))) - 0.5 * step;

    for (float x = min_x; x <= max_x; x += step) {
      for (float y = min_y; y <= max_y; y += step) {
        if (inside_triangle(x0, y0, x1, y1, x2, y2, x, y)) {
          Vector3D params = barycentric(x0, y0, x1, y1, x2, y2, x, y);
          int i = round((x - 0.5 * step) / step);
          int j = round((y - 0.5 * step) / step);
          Color color = params[0] * c0 + params[1] * c1 + params[2] * c2;
          sample_buffer[j * width * sqrt(sample_rate) + i] = color;
        }
      }
    }
  }



  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle

    const auto to_uv = [x0, y0, u0, v0,
                  x1, y1, u1, v1,
                  x2, y2, u2, v2](float x, float y) {
      Vector2D uv0 = Vector2D(u0, v0);
      Vector2D uv1 = Vector2D(u1, v1);
      Vector2D uv2 = Vector2D(u2, v2);
      Vector3D xy_barycentric = barycentric(x0, y0, x1, y1, x2, y2, x, y);
      return uv0 * xy_barycentric[0] + uv1 * xy_barycentric[1] + uv2 * xy_barycentric[2];
    };

    float step = 1.0 / sqrt(sample_rate);
    float min_x = floor(min(x0, min(x1, x2))) + 0.5 * step;
    float max_x = ceil(max(x0, max(x1, x2))) - 0.5 * step;
    float min_y = floor(min(y0, min(y1, y2))) + 0.5 * step;
    float max_y = ceil(max(y0, max(y1, y2))) - 0.5 * step;

    for (float x = min_x; x <= max_x; x += step) {
      for (float y = min_y; y <= max_y; y += step) {
        if (inside_triangle(x0, y0, x1, y1, x2, y2, x, y)) {
          Vector3D params = barycentric(x0, y0, x1, y1, x2, y2, x, y);
          int i = (int) (x / step);
          int j = (int) (y / step);
          SampleParams sp = SampleParams{to_uv(x, y),
                                          to_uv(x + 1, y),
                                          to_uv(x, y + 1),
                                          psm, lsm};
          Color color = tex.sample(sp);
          sample_buffer[j * width * sqrt(sample_rate) + i] = color;
        }
      }
    }
  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support
    this->sample_rate = rate;

    this->sample_buffer.resize(width * height * rate, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support
    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;


    this->sample_buffer.resize(width * height * this->sample_rate, Color::White);
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
        Color avg = Color(0, 0, 0);
        int super_x = x * sqrt(sample_rate);
        int super_y = y * sqrt(sample_rate);
        for (int i = super_x; i < super_x + sqrt(sample_rate); i++) {
          for (int j = super_y; j < super_y + sqrt(sample_rate); j++) {
            avg += (1.0 / sample_rate) * sample_buffer[j * width * sqrt(sample_rate) + i];
          }
        }
        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&avg.r)[k] * 255;
        }
      }
    }

  }

  Rasterizer::~Rasterizer() = default;


}// CGL
