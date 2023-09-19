# %%
import re
# %%
class GetSeqMarker:
    def __init__(self, parentheses_seq):
        self.parentheses_seq = parentheses_seq

    @property
    def getStMKs_premiRNA(self):
        marker = 0
        lstMKs_premiRNA = []
        cache = []
        for parenthese in self.parentheses_seq:
            if parenthese == "(":
                marker+=1
                lstMKs_premiRNA.append(marker)
                cache.append(marker)
            elif parenthese == ")":
                opposite_marker = cache.pop()
                lstMKs_premiRNA.append(opposite_marker)
            else:
                lstMKs_premiRNA.append(0)
        self.lstMKs_premiRNA = lstMKs_premiRNA

        list_site_markers = [[site, marker] for site, marker in enumerate(self.lstMKs_premiRNA)]
        self.list_site_markers = list_site_markers

    def getStMKs_miRNA3p(self,three_start_site, three_end_site):
        list_markers_three_miRNA= self.lstMKs_premiRNA[three_start_site:(three_end_site + 1)]
        self.list_markers_three_miRNA= list_markers_three_miRNA
        list_site_markers_three_miRNA = self.list_site_markers[three_start_site:(three_end_site + 1)]
        self.list_site_markers_three_miRNA= list_site_markers_three_miRNA

    def getStMKs_miRNA5p(self, five_start_site, five_end_site):
        list_markers_five_miRNA = self.lstMKs_premiRNA[five_start_site:(five_end_site + 1)]
        self.list_markers_five_miRNA = list_markers_five_miRNA
        list_site_markers_five_miRNA = self.list_site_markers[five_start_site:(five_end_site + 1)]
        self.list_site_markers_five_miRNA = list_site_markers_five_miRNA
    
    @property
    def rm0StMKs(self):
        try:
            lstMKsRM0_premiRNA = [marker for marker in self.lstMKs_premiRNA if marker!= 0]
            self.lstMKsRM0_premiRNA = lstMKsRM0_premiRNA
            
            lstStMKsRM0_premiRNA = [site_marker for site_marker in self.list_site_markers if site_marker[-1] != 0]
            self.lstStMKsRM0_premiRNA = lstStMKsRM0_premiRNA
        except:
            print("請先執行 .get_premiRNA_site_markers")
        try:
            lstMKsRM0_miRNAp3 = [marker for marker in self.list_markers_three_miRNA if marker != 0]
            lstMKsRM0_miRNAp3.reverse()
            self.lstMKsRM0_miRNAp3 = lstMKsRM0_miRNAp3
            
            lstStMKsRM0_miRNAp3 = [site_marker for site_marker in self.list_site_markers_three_miRNA if site_marker[-1] != 0]
            lstStMKsRM0_miRNAp3.reverse()
            self.lstStMKsRM0_miRNAp3 = lstStMKsRM0_miRNAp3
        except:
            print("請先執行 .get_three_miRNA")
        try:
            lstMKsRM0_miRNAp5 = [marker for marker in self.list_markers_five_miRNA if marker != 0]
            self.lstMKsRM0_miRNAp5 = lstMKsRM0_miRNAp5
            
            lstStMKsRM0_miRNAp5 = [site_marker for site_marker in self.list_site_markers_five_miRNA if site_marker[-1] != 0]
            self.lstStMKsRM0_miRNAp5 = lstStMKsRM0_miRNAp5
        except:
            print("請先執行 .get_five_miRNA")


# %%
class FindComparedSeq():
    def __init__(self, getMKs):
        self.lstMKs_premiRNA = getMKs.lstMKs_premiRNA
        self.lstStMKsRM0_premiRNA = getMKs.lstStMKsRM0_premiRNA
        self.lstStMKsRM0_miRNA = getMKs.lstStMKsRM0_miRNAp3
        #print('lstStMKsRM0_miRNA:',self.lstStMKsRM0_miRNA)
        self.startStRM0 = getMKs.lstStMKsRM0_miRNAp3[0][0]
        self.endStRM0 = getMKs.lstStMKsRM0_miRNAp3[-1][0]
        self.startMKsRM0 = getMKs.lstStMKsRM0_miRNAp3[0][1]
        self.endMKsRM0 = getMKs.lstStMKsRM0_miRNAp3[-1][1]
        '''
        print('self.startStRM0:',self.startStRM0)
        print('self.endStRM0:',self.endStRM0)
        print('self.startMKsRM0:',self.startMKsRM0)
        print('self.endMKsRM0:',self.endMKsRM0)
        self.list_markers_end = getMKs.lstMKsRM0_miRNAp3
        '''
        self.p53 = '3p'
        

    @property
    def get_index_c(self):
        final_index = 0
        first_index = 0
        for index, site_markers_rm0 in enumerate(self.lstStMKsRM0_premiRNA):
            site = site_markers_rm0[0]
            marker = site_markers_rm0[1]
            if marker == self.startMKsRM0 and site != self.startStRM0:
                first_index = index
                self.first_index = first_index
            elif marker == self.endMKsRM0 and site != self.endStRM0:
                final_index = index
                self.final_index = final_index
        '''
        print('first_index:', first_index)
        print('final_index:', final_index)
        '''

    @property
    def get_start_site_c(self):
        # if self.startMKsRM0 == self.lstStMKsRM0_premiRNA[0][1] and self.startStRM0 in getMKs.lstMKsRM0_miRNAp3:
        # print('self:', self.lstStMKsRM0_premiRNA)
        if self.p53 == '3p':
            if self.startMKsRM0 == self.lstStMKsRM0_premiRNA[0][1]:
                # print(self.startMKsRM0)
                # print(self.lstStMKsRM0_premiRNA[0][1])
                start_site_c = 0
            else:
                prev_index = self.first_index - 1
                # print('prev_index:',prev_index)
                start_site_c = self.lstStMKsRM0_premiRNA[prev_index][0] + 1
            self.start_site_c = start_site_c
        elif self.p53 == '5p':
            '''
            print('w:',self.startMKsRM0)
            print('x:',self.lstStMKsRM0_premiRNA[-1][1])
            '''
            if self.startMKsRM0 == self.lstStMKsRM0_premiRNA[-1][1] and self.startMKsRM0 in self.lstStMKsRM0_miRNA:
                # print(self.startMKsRM0)
                # print(self.lstStMKsRM0_premiRNA[0][1])
                start_site_c = 0
            else:
                prev_index = self.first_index - 1
                # print('prev_index:',prev_index)
                start_site_c = self.lstStMKsRM0_premiRNA[prev_index][0] + 1
            self.start_site_c = start_site_c
    @property
    def get_end_site_c(self):
        # if self.endMKsRM0 == self.lstStMKsRM0_premiRNA[-1][1] and self.endStRM0 in getMKs.lstMKsRM0_miRNAp5:
        # print('self:', self.lstStMKsRM0_premiRNA)
        '''
        print('first_index', self.first_index)
        print('final_index', self.final_index)
        print(self.endMKsRM0)
        print(self.startMKsRM0)
        print(self.lstStMKsRM0_premiRNA)
        print(self.lstStMKsRM0_premiRNA[-1][1])
        print(self.p53)
        '''
        if self.p53 == '3p':
            if self.endMKsRM0 == self.lstStMKsRM0_premiRNA[-1][1]:
                end_site_c =  len(self.lstMKs_premiRNA)
            else:
                next_index = self.final_index + 1
                # print('next_index:',next_index)
                end_site_c = self.lstStMKsRM0_premiRNA[next_index][0] - 1
            self.end_site_c = end_site_c
            # print(self.end_site_c, self.start_site_c)
        if self.p53 == '5p':
            if self.startMKsRM0 == self.lstStMKsRM0_premiRNA[-1][1]:
                end_site_c =  len(self.lstMKs_premiRNA)
            else:
                next_index = self.final_index + 1
                # print('next_index:',next_index)
                end_site_c = self.lstStMKsRM0_premiRNA[next_index][0] - 1
            self.end_site_c = end_site_c
            # print(self.start_site_c, self.end_site_c)

    @property
    def get_list_compared_markers(self):
        compared_markers = self.lstMKs_premiRNA[self.start_site_c:(self.end_site_c+1)]
        self.compared_markers = compared_markers
        # print('start_site:', self.start_site_c, 'end_site+1:', self.end_site_c+1)
        #compared_markers.reverse()
        # return compared_markers

# %%
FindComparedSeq.mro()
# %%
class ThreeEnd_miRNA(FindComparedSeq):
    pass

# %%  
class FiveEnd_miRNA(FindComparedSeq):
    def __init__(self, getMKs):
        super().__init__(getMKs)
        self.lstStMKsRM0_miRNA = getMKs.lstStMKsRM0_miRNAp5
        # print('lstStMKsRM0_miRNA:',self.lstStMKsRM0_miRNA)
        self.startStRM0 = getMKs.lstStMKsRM0_miRNAp5[0][0]
        self.endStRM0 = getMKs.lstStMKsRM0_miRNAp5[-1][0]
        self.startMKsRM0 = getMKs.lstStMKsRM0_miRNAp5[0][1]
        self.endMKsRM0 = getMKs.lstStMKsRM0_miRNAp5[-1][1]
        '''
        print('self.startStRM0:',self.startStRM0)
        print('self.endStRM0:',self.endStRM0)
        print('self.startMKsRM0:',self.startMKsRM0)
        print('self.endMKsRM0:',self.endMKsRM0)
        '''
        self.list_markers_end = getMKs.lstMKsRM0_miRNAp5
        self.p53 = '5p'

    @property
    def get_index_c(self):
        first_index = 0
        final_index = 0
        for index, site_markers_rm0 in enumerate(self.lstStMKsRM0_premiRNA):
            site = site_markers_rm0[0]
            marker = site_markers_rm0[1]
            if marker == self.endMKsRM0 and site != self.endStRM0:
                first_index = index
                self.first_index = first_index
            elif marker == self.startMKsRM0 and site != self.startStRM0:
                final_index = index
                #print('elif marker == self.startMKsRM0 and site != self.startStRM0:')
                #print(f'elif {marker} == {self.startMKsRM0} and {site} != {self.startStRM0}:')
                #print('final_index', final_index)
                self.final_index = final_index
### 判斷字串是否連續 (FUNCTION)
def cotinous_seq(list_markers):
    raw_markers = [num for num in list_markers if num != 0]
    if len(raw_markers) == 0:
        return None

    sign = None
    plus_sign = +1
    minus_sign = -1

    result_seq_markers = [[]]
    for num in raw_markers:
        if len(result_seq_markers[-1]) == 0:
            result_seq_markers[-1].append(num)
        elif len(result_seq_markers[-1]) == 1:
            if num == result_seq_markers[-1][-1] + 1:
                sign = plus_sign
            elif num == result_seq_markers[-1][-1] - 1:
                sign = minus_sign
            else:
                result_seq_markers.append([])
            result_seq_markers[-1].append(num)
        else:
            if num == result_seq_markers[-1][-1] + sign:
                pass
            else:
                result_seq_markers.append([])
            result_seq_markers[-1].append(num)
    return result_seq_markers
    