import numpy as np


class match_1d():
    def __init__(self,set1,set2) -> None:
        self.distance_error = 0.4
        self.n1 = len(set1)
        self.n2 = len(set2)
        self.cor  = np.zeros((self.n1,self.n2),dtype=np.double)
        self.shift_array = np.zeros(self.n1*self.n2,dtype=np.double)

        self.shift_cor(set1,set2)
        self.get_best_shift()

    
    def shift_cor(self,set1,set2):
        for i in range(self.n1):
            for j in range(self.n2):
                dis = set2[j] - set1[i]
                self.cor[i,j] = dis
                self.shift_array[i*self.n2+j] = dis
    
    def get_best_shift(self):
        t_array = np.sort(self.shift_array)
        t_ids = np.argsort(self.shift_array)
        set1_id = t_ids//self.n2
        set2_id = t_ids%self.n2
        # 双指针
        fp = -1
        ep = 0
        #平移量和对应的个数
        t_message = []
        t_number = 0
        while ep < t_array.shape[0]:
            fp += 1
            if t_array[ep] - t_array[fp] > 2 * self.distance_error:
                continue
            while ep < t_array.shape[0] and t_array[ep] -  t_array[fp] <= 2 * self.distance_error:
                ep += 1
            #fp to ep-1
            ids1 = len(set(set1_id[fp:ep]))
            ids2 = len(set(set2_id[fp:ep]))
            now_number = min(ids1,ids2)
            print((now_number,fp,ep,t_array[fp]))
            if now_number > t_number:
                t_message = (now_number,fp,ep)
                t_number = now_number

        print(t_message)
                

if __name__ == "__main__":
    old = [1,2,3,4,5]
    new = [6,7,8,9,11,12]
    match_1d(old,new)