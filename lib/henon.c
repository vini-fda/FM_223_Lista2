#define THRESHOLD 1.0e3

float iterated_henon_map(float a, float b, float x, float y, int n) {
    float x_ = x;
    for (int i = 0; i < n; i++) {
        x_ = a - x * x + b * y;
        y = x;
        x = x_;
        // return 1.0 if the orbit escapes
        if (x > THRESHOLD || x < -THRESHOLD || y > THRESHOLD || y < -THRESHOLD) {
            return 1.0;
        }
    }
    return 0.0;
}