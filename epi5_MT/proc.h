#pragma once
#include <string>
void proc(const std::string& init_data_path, bool force_cornif,bool use_last=false,const std::string& init_uvp_data="",const std::string& init_w_data="",const std::string& init_ATP_data="",
          const std::string& init_ext_stim_data="");
