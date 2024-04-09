#include "lowpass_filter.h"

LowPassFilter::LowPassFilter(float time_constant)
    : Tf(time_constant)
    , y_prev(0.0f)
{
    timestamp_prev = _micros();
}


float LowPassFilter::operator() (float x)
{
    unsigned long timestamp = _micros();
    float dt = (timestamp - timestamp_prev)*1e-6f;

    if (dt < 0.0f ) dt = 1e-3f;
    else if(dt > 0.3f) {
        y_prev = x;
        timestamp_prev = timestamp;
        return x;
    }
    timestamp_prev = timestamp;
    return calc_filter(x,dt);
}

float LowPassFilter::operator() (float x, float Ts)
{
    return calc_filter(x,Ts);
}

inline float LowPassFilter::calc_filter(float x, float dt)
{
    float alpha = Tf/(Tf + dt);
    float y = alpha*y_prev + (1.0f - alpha)*x;
    y_prev = y;
    return y;
}
