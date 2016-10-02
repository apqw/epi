#include "color.h"
#include "../global.h"
#include "VisParams.h"
#include <vector>
#include "vis_global.h"
#include "../calc/cpu2/CellManager.h"
#include <algorithm>
static const std::vector<float> c_pk = { 1.0f,0.6f,0.8f };
static const std::vector<float> c_white = { 1.0f,1.0f,1.0f };
static const std::vector<float> c_pkwhite = { 1.0f,0.8f,1.0f };
static const std::vector<float> c_memb = { 0.5f,0.5f,0.5f };
static const std::vector<float> c_fix_malig = { 0.8f,0.0f,0.1f };
static const std::vector<float> c_fix_normal = { 0.0f,0.8f,0.1f };


static const constexpr double UMIN = 0.1f, UMAX = 0.65f;
static const constexpr double FMIN = 0.0f, FMAX = 1.0f;
template <typename T>
static bool is_in_bounds(const T& value, const T& low, const T& high) {
    return !(value < low) && (value < high); //test [low,high)
}

static double clamp(const double value, const double min, const double max) {
    return std::max(std::min(value, max), min);
}
static double lerp_alpha(const double value, const double min, const double max) {
    return (clamp(value, min, max) - min) / (max - min);
}
static std::vector<float> c_color_M(float c, unsigned int fix_orig) {
    const bool is_malig = fix_orig < pm->MALIG_NUM;
    const std::vector<float> malig_color0 = { 1.0f,0.5f,0.5f };
    const std::vector<float> malig_color1 = { 1.0f,0.0f,0.0f };

    const std::vector<float> normal_color_low = { 1.0f,0.6f,0.8f };
    const std::vector<float> normal_color1 = { 0.5f,c,0.5f };
    const std::vector<float> normal_color2 = { 0.0f,1.0f,2.0f - 4.0f*c };
    const std::vector<float> normal_color3 = { 4.0f*c - 2.0f,1.0f,0.0f };
    const std::vector<float> normal_color4 = { 1.0f,4.0f - 4.0f*c,0.0f };
    const std::vector<float> normal_color_high = { 1.0f,0.0f,0.0f };

    const std::vector<const std::vector<float>*> normal_color_set_in_order
        = { &normal_color1,&normal_color2,&normal_color3,&normal_color4 };
    const int c_num = (int)normal_color_set_in_order.size();
    if (is_malig) {
        if (c < 0.0) {
            return malig_color0;
        }
        else {
            return malig_color1;
        }
    }
    else {
        if (c < 0.0) {
            return normal_color_low;
        }
        if (!(c < 1.0)) {
            return normal_color_high;
        }

        for (size_t i = 0; i < c_num; i++) {
            if (is_in_bounds(c, i / (float)c_num, (i + 1) / (float)c_num)) {
                return *normal_color_set_in_order[i];
            }
        }

        throw std::logic_error("An uncovered range exists.[1]");

    }
}


static std::vector<float> c_color_default(float c) {

    const std::vector<float> normal_color_low = { 0.0f,0.0f,1.0f };
    const std::vector<float> normal_color1 = { 0.0f,4.0f*c,1.0f };
    const std::vector<float> normal_color2 = { 0.0f,1.0f,2.0f - 4.0f*c };
    const std::vector<float> normal_color3 = { 4.0f*c - 2.0f,1.0f,0.0f };
    const std::vector<float> normal_color4 = { 1.0f,4.0f - 4.0f*c,0.0f };
    const std::vector<float> normal_color_high = { 1.0f,0.0f,0.0f };

    const std::vector<const std::vector<float>*> normal_color_set_in_order
        = { &normal_color1,&normal_color2,&normal_color3,&normal_color4 };
    const int c_num = (int)normal_color_set_in_order.size();

    if (c < 0.0) {
        return normal_color_low;
    }
    if (!(c < 1.0)) {
        return normal_color_high;
    }

    for (size_t i = 0; i < c_num; i++) {
        if (is_in_bounds(c, i / (float)c_num, (i + 1) / (float)c_num)) {
            return *normal_color_set_in_order[i];
        }
    }

    throw std::logic_error("An uncovered range exists.[2]");


}

static std::vector<float> c_color_musume(unsigned int div, unsigned int tb)
{
    const bool is_malig = tb < pm->MALIG_NUM;

    static const std::vector<float> malig_musume = { 1.0f,0.7f,0.7f };
    static const std::vector<float> normal_musume = { 0.7f,1.0f,0.7f };

    static const std::vector<float> div_add = { 0.2f,0.2f,0.2f };
    std::vector<float> tmp;
    if (is_malig) {
        tmp = malig_musume;
    }
    else {
        tmp = normal_musume;
    }

    if (div == 0) {
        std::transform(tmp.begin(), tmp.end(), div_add.begin(), tmp.begin(), std::plus<float>());
    }

    return tmp;

}

static std::vector<float> c_color_FIX(unsigned int fix_orig) {
    const bool is_malig = fix_orig < pm->MALIG_NUM;
    if (is_malig) {
        return c_fix_malig;
    }
    else {
        return c_fix_normal;
    }
}
static std::vector<float> get_MEMB_color(const CellTempStruct&cts) {
    switch (vp.mode)
    {
    case VisParams::MEDICAL:
        return c_pkwhite;
        break;
    default:
        return c_memb;
        break;
    }
}

static std::vector<float> get_FIX_color(const CellTempStruct&cts) {
    switch (vp.mode)
    {
    case VisParams::MEDICAL:
        return c_color_M(-1.0f, cts.stem_orig_id);
        break;
    default:
        return c_color_FIX(cts.stem_orig_id);
        break;
    }
}

static std::vector<float> ca2p_color(const CellTempStruct&cts) {
    return c_color_default((float)lerp_alpha(cts.ca2p_avg, UMIN, UMAX));
}

static std::vector<float> ex_fat_color(const CellTempStruct&cts) {
    return c_color_default((float)lerp_alpha(cts.ex_fat, FMIN, FMAX));
}

static std::vector<float> in_fat_color(const CellTempStruct&cts) {
    return c_color_default((float)lerp_alpha(cts.fat, FMIN, FMAX));
}

static std::vector<float> color_modewise(const CellTempStruct&cts) {
    switch (vp.mode)
    {
    case VisParams::MEDICAL:
        return c_color_M(0.122f, cts.stem_orig_id);
    case VisParams::CA2P:
        return ca2p_color(cts);
    case VisParams::EX_FAT:
        return ex_fat_color(cts);
    case VisParams::IN_FAT:
        return in_fat_color(cts);
    default:
        throw std::runtime_error("Unknown display mode:"_s + std::to_string(vp.mode));
    }
}
static std::vector<float> get_MUSUME_color(const CellTempStruct&cts) {

    switch (vp.disp_musume) {
    case 1:
        return color_modewise(cts);
    case 2:
        return c_color_musume(cts.div_times, cts.stem_orig_id);
    default:
        throw std::runtime_error("Unknown MUSUME display mode:"_s+std::to_string(vp.disp_musume));
    }
}

static std::vector<float> get_ALIVE_color(const CellTempStruct&cts) {
    return color_modewise(cts);
}

static std::vector<float> get_DEAD_color(const CellTempStruct&cts) {
    switch (vp.mode) {
    case VisParams::MEDICAL:
        return c_pk;
    case VisParams::CA2P:
        return c_white;
    case VisParams::EX_FAT:
        return ex_fat_color(cts);
    case VisParams::IN_FAT:
        return in_fat_color(cts);
    default:
        throw std::runtime_error("Unknown display mode:"_s + std::to_string(vp.mode));
    }
}

std::vector<float> get_color_rgb(const CellTempStruct&cts) {
    switch (cts.state) {
    case MEMB:
        return get_MEMB_color(cts);
    case FIX:
        return get_FIX_color(cts);
    case MUSUME:
        return get_MUSUME_color(cts);
    case ALIVE:
        return get_ALIVE_color(cts);
    case DEAD:
        return get_DEAD_color(cts);
    default:
        return{ 1.0f,1.0f,1.0f };
        break;
    }
}