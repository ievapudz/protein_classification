#ifndef _CONSTANTS_
#define _CONSTANTS_

class Constants{
    private:
        const static double gap_open_penalty_ = 4.0;
        const static double gap_ext_penalty_ = 1.0;
    public:
        static double gapOpenPenalty(){
            return gap_open_penalty_;
        }
        static double gapExtPenalty(){
            return gap_ext_penalty_;
        }
};

#endif
