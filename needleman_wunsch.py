
class Coordinates:
    def __str__(self):
        return '({},{})'.format(self.col, self.row)

class Data:
    
    def __init__(self):
        self._current = Coordinates()
        self._prev = Coordinates()

    def __str__(self):
        return '{}'.format(self._value)

class Dictionary:

    def __init__(self):
        # Create Dictionary
        file_path = 'BLOSUM62.csv'
        labels = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X','*']
        self._dictionary = {}
        try:
            file = open(file_path)
            for i, line in enumerate(file):
                self._dictionary[labels[i]] = {}
                for j, attr in enumerate(line.strip().split(',')):
                    self._dictionary[labels[i]][labels[j]] = int(attr)
        finally:
            file.close()


class Alignment(Dictionary):
    
    def __init__(self):
        super().__init__()

    def matching(self, sequence_1, sequence_2):
        self.__sequence_1 = sequence_1
        self.__sequence_2 = sequence_2
        self.__rows = len(sequence_1)
        self.__cols = len(sequence_2)
        self.__init_matrix()
        self.__find_max_value()
        self.__show_matrix()
        self.__new_sequence()

    def __show_matrix(self):
        # Print score
        for i in range(self.__rows + 1):
            for j in range(self.__cols + 1):
                print(self.__matrix[i][j], end='\n' if j == self.__cols else '\t')
        print()
        # Print Prev coor
        for i in range(self.__rows + 1):
            for j in range(self.__cols + 1):
                print(self.__matrix[i][j]._prev if not(i==0 or j==0) else None, end='\n' if j == self.__cols else '\t')

    def __new_sequence(self):
        new_sequence_1 = ''
        new_sequence_2 = ''

        print('score :', self.__maxi._value)
        current = self.__maxi._current
        prev = self.__maxi._prev
        
        # Add post indels sequence_1
        for i in range(self.__cols - current.col):
            new_sequence_1 = '-' + new_sequence_1
            new_sequence_2 = self.__sequence_2[self.__cols - i - 1] + new_sequence_2
	    # Add post indels sequence_2
        for i in range(self.__rows - current.row):
            new_sequence_1 = self.__sequence_1[self.__rows - i - 1] + new_sequence_1
            new_sequence_2 = '-' + new_sequence_2        

        # Traceback
        while self.__maxi._current is not None:
            direction = self.__direction()
            current = self.__maxi._current
            prev = self.__maxi._prev
            if direction == 'row':
                distance = current.col - prev.col
                for i in range(1, distance + 1):
                    new_sequence_1 = '-' + new_sequence_1
                    new_sequence_2 = self.__sequence_2[current.col - i] + new_sequence_2
            elif direction == 'column':
                distance = current.row - prev.row
                for i in range(1, distance + 1):
                    new_sequence_1 = self.__sequence_1[current.row - i] + new_sequence_1
                    new_sequence_2 = '-' + new_sequence_2
            else:
                new_sequence_1 = self.__sequence_1[prev.row] + new_sequence_1
                new_sequence_2 = self.__sequence_2[prev.col] + new_sequence_2
            self.__maxi = self.__matrix[prev.row][prev.col]

	    # Add pre indels sequence_1
        for i in range(1, current.col):
            new_sequence_1 = '-' + new_sequence_1
            new_sequence_2 = self.__sequence_2[current.col - i - 1] + new_sequence_2
        # Add pre indels sequence_2
        for i in range(1, current.row):
            new_sequence_1 = self.__sequence_1[current.row - i - 1] + new_sequence_1
            new_sequence_2 = '-' + new_sequence_2

        print(new_sequence_1)
        print(new_sequence_2)

    def __direction(self):
        if self.__maxi._current.row == self.__maxi._prev.row:
            return 'row'
        elif self.__maxi._current.col == self.__maxi._prev.col:
            return 'column'
        else:
            return 'diagonal'

    def __find_max_value(self):
        rows = len(self.__sequence_1)
        cols = len(self.__sequence_2)

        self.__maxi = self.__matrix[rows][cols]

        # Max in a row
        for i in range(1, cols):
            temp = self.__matrix[rows][i]._value + (-6 * (cols - i))
            if self.__maxi._value < temp:
                self.__maxi = self.__matrix[rows][i]
                self.__maxi._value = temp

        # Max in a column
        for i in range(1, rows):
            temp = self.__matrix[i][cols]._value + (-6 * (rows - i))
            if self.__maxi._value < temp:
                self.__maxi = self.__matrix[i][cols]
                self.__maxi._value = temp

    def __init_matrix(self):
        # initial self.__matrix
        # Needleman-Wunsch Algorithm
        self.__matrix = []
        for i in range(len(self.__sequence_1) + 1):
            self.__matrix.append([])
            for j in range(len(self.__sequence_2) + 1):
                data = Data()
                if i == 0 or j == 0:
                    if i == j:
                        data._value = 0
                    else:
                        data._value = -6 * (j if i == 0 else i)
                    data._current = None
                else:
                    cross = self.__matrix[i - 1][j - 1]._value + self._dictionary[self.__sequence_1[i - 1]][self.__sequence_2[j - 1]]
                    # Max value in a row
                    row, path_col = self.__max_value([self.__matrix[i][col]._value for col in range(j)])
                    # Max value in a column
                    col, path_row = self.__max_value([self.__matrix[row][j]._value for row in range(i)])

                    if cross >= row and cross >= col:
                        data._value = cross
                        data._prev.row = i - 1
                        data._prev.col = j - 1
                    elif row >= col:
                        data._value = row
                        data._prev.row = i
                        data._prev.col = path_col
                    else:
                        data._value = col
                        data._prev.row = path_row
                        data._prev.col = j

                    data._current.row = i
                    data._current.col = j
                self.__matrix[i].append(data)

    def __max_value(self, values):
        n = len(values)
        maxi = -1000000000
        path = 0
        for i, value in enumerate(values):
            max_value = value + (-6 * (n - i))
            if maxi < max_value:
                maxi = max_value
                path = i
        return maxi, path

def main():
    algin = Alignment()
    sequence_1 = 'QALVAYA'
    sequence_2 = 'NALWVAYMA'
    print(sequence_1)
    print(sequence_2)
    # sequence_1 = input()
    # sequence_2 = input()
    algin.matching(sequence_1, sequence_2)

if __name__ == '__main__':
    main()
