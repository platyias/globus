//
//  z_score_map.h
//
//  Created by Jeewoen Shin on 4/14/17.
//

#include <cmath>

#include <map>
using std::map;

#include <string>
using std::string;

#ifndef z_score_map_h
#define z_score_map_h

void map_gcZ_Pnb (map<int, double> &gcZ_Pnb){
    // score -> probability -> log(probability)
    gcZ_Pnb[ 0] = log10(0.0029);   // 0-0.5
    gcZ_Pnb[ 1] = log10(0.0029);   // 0.5-1
    gcZ_Pnb[ 2] = log10(0.0036);   // 1-1.5
    gcZ_Pnb[ 3] = log10(0.0036);   // 1.5-2
    gcZ_Pnb[ 4] = log10(0.0050);   // 2-2.5
    gcZ_Pnb[ 5] = log10(0.0050);   // 2.5-3
    gcZ_Pnb[ 6] = log10(0.0053);   // 3-3.5
    gcZ_Pnb[ 7] = log10(0.0053);   // 3.5-4**************0.0153
    gcZ_Pnb[ 8] = log10(0.0087);   // 4-4.5
    gcZ_Pnb[ 9] = log10(0.0087);   // 4.5-5
    gcZ_Pnb[10] = log10(0.0089);   // 5-5.5
    gcZ_Pnb[11] = log10(0.0089);   // 5.5-6 *0.2407
    gcZ_Pnb[12] = log10(0.0078);   // 6-6.5 *0.2467
    gcZ_Pnb[13] = log10(0.0078);   // 6.5-7 *0.2167
    gcZ_Pnb[14] = log10(0.0312);   // 7-7.5 *0.2437
    gcZ_Pnb[15] = log10(0.0312);   // 7.5-8
    gcZ_Pnb[16] = log10(0.0614);   // 8-8.5
    gcZ_Pnb[17] = log10(0.0614);   // 8.5-9
    gcZ_Pnb[18] = log10(0.1129);   // 9-9.5
    gcZ_Pnb[19] = log10(0.1129);   // 9.5-10
    gcZ_Pnb[20] = log10(0.1429);   // 10-10.5
    gcZ_Pnb[21] = log10(0.1429);   // 10.5-11
    gcZ_Pnb[22] = log10(0.1429);   // 11-11.5
    gcZ_Pnb[23] = log10(0.1429);   // 11.5-12
    gcZ_Pnb[24] = log10(0.1970);   // 12-12.5
    gcZ_Pnb[25] = log10(0.1970);   // 12.5-13
    gcZ_Pnb[26] = log10(0.1970);   // 13-13.5
    gcZ_Pnb[27] = log10(0.1970);   // 13.5-14
    gcZ_Pnb[28] = log10(0.1786);   // 14-14.5
    gcZ_Pnb[29] = log10(0.1786);   // 14.5-15
    gcZ_Pnb[30] = log10(0.1786);   // 15-15.5
    gcZ_Pnb[31] = log10(0.1786);   // 15.5-16
    gcZ_Pnb[32] = log10(0.2750);   // 16-16.5
    gcZ_Pnb[33] = log10(0.2750);   // 16.5-17
    gcZ_Pnb[34] = log10(0.2750);   // 17-17.5
    gcZ_Pnb[35] = log10(0.2750);   // 17.5-18
    gcZ_Pnb[36] = log10(0.3043);   // 18-18.5
    gcZ_Pnb[37] = log10(0.3043);   // 18.5-19
    gcZ_Pnb[38] = log10(0.3043);   // 19-19.5
    gcZ_Pnb[39] = log10(0.3043);   // 19.5-20
    gcZ_Pnb[40] = log10(0.6418);   // >=20

}

void map_phZ_Pnb (map<int, double> &phZ_Pnb){
    // score -> probability -> log(probability)
    phZ_Pnb[ 0] = log10(0.0028);   // 0-0.5
    phZ_Pnb[ 1] = log10(0.0026);   // 0.5-1
    phZ_Pnb[ 2] = log10(0.0038);   // 1-1.5
    phZ_Pnb[ 3] = log10(0.0041);   // 1.5-2
    phZ_Pnb[ 4] = log10(0.0057);   // 2-2.5
    phZ_Pnb[ 5] = log10(0.0056);   // 2.5-3
    phZ_Pnb[ 6] = log10(0.0087);   // 3-3.5
    phZ_Pnb[ 7] = log10(0.0143);   // 3.5-4
    phZ_Pnb[ 8] = log10(0.0149);   // 4-4.5
    phZ_Pnb[ 9] = log10(0.0574);   // 4.5-5
    phZ_Pnb[10] = log10(0.0930);   // 5-5.5
    phZ_Pnb[11] = log10(0.2407);   // 5.5-6 *0.2407
    phZ_Pnb[12] = log10(0.2437);   // 6-6.5 *0.2467
    phZ_Pnb[13] = log10(0.2437);   // 6.5-7 *0.2167
    phZ_Pnb[14] = log10(0.2467);   // 7-7.5 *0.2437
    phZ_Pnb[15] = log10(0.2467);   // 7.5-8
    phZ_Pnb[16] = log10(0.2467);   // 8-8.5
    phZ_Pnb[17] = log10(0.2467);   // 8.5-9
    phZ_Pnb[18] = log10(0.2467);   // 9-9.5
    phZ_Pnb[19] = log10(0.2467);   // 9.5-10
    phZ_Pnb[20] = log10(0.2467);   // 10-10.5
    phZ_Pnb[21] = log10(0.2467);   // 10.5-11
    phZ_Pnb[22] = log10(0.2467);   // 11-11.5
    phZ_Pnb[23] = log10(0.2467);   // 11.5-12
    phZ_Pnb[24] = log10(0.2467);   // 12-12.5
    phZ_Pnb[25] = log10(0.2467);   // 12.5-13
    phZ_Pnb[26] = log10(0.2467);   // 13-13.5
    phZ_Pnb[27] = log10(0.2467);   // 13.5-14
    phZ_Pnb[28] = log10(0.2467);   // 14-14.5
    phZ_Pnb[29] = log10(0.2467);   // 14.5-15
    phZ_Pnb[30] = log10(0.2467);   // 15-15.5
    phZ_Pnb[31] = log10(0.2467);   // 15.5-16
    phZ_Pnb[32] = log10(0.2467);   // 16-16.5
    phZ_Pnb[33] = log10(0.2467);   // 16.5-17
    phZ_Pnb[34] = log10(0.2467);   // 17-17.5
    phZ_Pnb[35] = log10(0.2467);   // 17.5-18
    phZ_Pnb[36] = log10(0.2467);   // 18-18.5
    phZ_Pnb[37] = log10(0.2467);   // 18.5-19
    phZ_Pnb[38] = log10(0.2467);   // 19-19.5
    phZ_Pnb[39] = log10(0.2467);   // 19.5-20
    phZ_Pnb[40] = log10(0.2467);   // >=20
    

}

void map_ceZ_Pnb (map<int, double> &ceZ_Pnb){
    // score -> probability -> log(probability)
    ceZ_Pnb[ 0] = log10(0.0024);   // 0-0.5
    ceZ_Pnb[ 1] = log10(0.0027);   // 0.5-1
    ceZ_Pnb[ 2] = log10(0.0028);  // 1-1.5
    ceZ_Pnb[ 3] = log10(0.0027);   // 1.5-2
    ceZ_Pnb[ 4] = log10(0.0041);   // 2-2.5
    ceZ_Pnb[ 5] = log10(0.0076);   // 2.5-3
    ceZ_Pnb[ 6] = log10(0.0154); // 3-3.5
    ceZ_Pnb[ 7] = log10(0.0275);  // 3.5-4
    ceZ_Pnb[ 8] = log10(0.0387); // 4-4.5
    ceZ_Pnb[ 9] = log10(0.0790);   // 4.5-5
    ceZ_Pnb[10] = log10(0.1455); // 5-5.5
    ceZ_Pnb[11] = log10(0.1824);   // 5.5-6 *0.2407
    ceZ_Pnb[12] = log10(0.2727); // 6-6.5 *0.2467
    ceZ_Pnb[13] = log10(0.2041);   // 6.5-7 *0.2167
    ceZ_Pnb[14] = log10(0.3889); // 7-7.5 *0.2437
    ceZ_Pnb[15] = log10(0.3889);   // 7.5-8
    ceZ_Pnb[16] = log10(0.3889); // 8-8.5
    ceZ_Pnb[17] = log10(0.3889);   // 8.5-9
    ceZ_Pnb[18] = log10(0.3889); // 9-9.5
    ceZ_Pnb[19] = log10(0.3889);   // 9.5-10
    ceZ_Pnb[20] = log10(0.3889);  // 10-10.5
    ceZ_Pnb[21] = log10(0.3889);   // 10.5-11
    ceZ_Pnb[22] = log10(0.3889); // 11-11.5
    ceZ_Pnb[23] = log10(0.3889);   // 11.5-12
    ceZ_Pnb[24] = log10(0.3889); // 12-12.5
    ceZ_Pnb[25] = log10(0.3889);   // 12.5-13
    ceZ_Pnb[26] = log10(0.3889); // 13-13.5
    ceZ_Pnb[27] = log10(0.3889);   // 13.5-14
    ceZ_Pnb[28] = log10(0.3889); // 14-14.5
    ceZ_Pnb[29] = log10(0.3889);   // 14.5-15
    ceZ_Pnb[30] = log10(0.3889);// 15-15.5
    ceZ_Pnb[31] = log10(0.3889);   // 15.5-16
    ceZ_Pnb[32] = log10(0.3889); // 16-16.5
    ceZ_Pnb[33] = log10(0.3889);   // 16.5-17
    ceZ_Pnb[34] = log10(0.3889); // 17-17.5
    ceZ_Pnb[35] = log10(0.3889);   // 17.5-18
    ceZ_Pnb[36] = log10(0.3889); // 18-18.5
    ceZ_Pnb[37] = log10(0.3889);   // 18.5-19
    ceZ_Pnb[38] = log10(0.3889); // 19-19.5
    ceZ_Pnb[39] = log10(0.3889);   // 19.5-20
    ceZ_Pnb[40] = log10(0.3889); // >=20
}

void map_hmZ_Ptr (map<int, double> &hmZ_Ptr){
    // score -> probability -> log(probability)
    hmZ_Ptr[0] = -4.0;             // 0-5%
    hmZ_Ptr[1] = -4.0;             // 5-10%
    hmZ_Ptr[2] = -4.0;             // 10-15%
    hmZ_Ptr[3] = -4.0;             // 15-20%
    hmZ_Ptr[4]  = log10(0.0091);  // 20-25% // -2.04096
    hmZ_Ptr[5]  = log10(0.0268);  // 25-30%
    hmZ_Ptr[6]  = log10(0.0554);  // 30-35%
    hmZ_Ptr[7]  = log10(0.1160);  // 35-40%
    hmZ_Ptr[8]  = log10(0.3345);  // 40-45%
    hmZ_Ptr[9]  = log10(0.3882);  // 45-50%
    hmZ_Ptr[10] = log10(0.5664);  // 50-55%
    hmZ_Ptr[11] = log10(0.6629);  // 55-60%
    hmZ_Ptr[12] = log10(0.8438);  // 60-65%
    hmZ_Ptr[13] = log10(0.8514);  // 65-70%
    hmZ_Ptr[14] = log10(0.9130);  // 70-75%
    hmZ_Ptr[15] = log10(0.9224);  // 75-80%
    hmZ_Ptr[16] = log10(0.9443);  // 80-85%
    hmZ_Ptr[17] = log10(0.9487);  // 85-90%
    hmZ_Ptr[18] = log10(0.9706);  // 90-95%
    hmZ_Ptr[19] = log10(0.9750);  // 95-100%
}

#endif
