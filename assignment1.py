"""
Assignment 1 for the course UE18MA251, Linear Algebra

Given a matrix A for the system of equations Ax=b
    ->Display, number of equations/rows, number of unknowns/columns, and rank of the matrix. (Use built in functions if required)
    ->Is the matrix singular (if it is a square matrix)?

    ->For a given b, will the system Ax=b have

        ->Unique solution.
        ->Infinitely many solutions
        ->No Solution.

    ->If it has a unique solution display the solution.
"""

class LinearSystem:
    def __init__(self):
        pass
        
    def get_det(self, A):
        if not self.is_square(A):
            print("Determinant of non-square matrices not possible")
        else:
            if len(A)==2:
                return A[0][0]*A[1][1]-A[0][1]*A[1][0]
            det = 0
            for c in range(len(A)):
                det += ((-1)**c)*A[0][c]*self.get_det(self.get_minor_matrix(A,0,c))
                
            return det

    def get_minor_matrix(self, A, r, c):
        return [row[:c] + row[c+1:] for row in (A[:r] + A[r+1:])]
    
    def get_zero_matrix(self, r, c):
        return [[0 for _ in range(c)] for _ in range(r)]
    
    def get_dim(self, M):
        return (len(M), len(M[0]))
                    
    def copy_matrix(self, M):
        dim = self.get_dim(M)
        rows = dim[0]
        cols = dim[1]
        copy = self.get_zero_matrix(rows, cols)
        
        for i in range(rows):
            for j in range(cols):
                copy[i][j] = M[i][j]
                
        return copy
    
    def is_square(self, M):
        dim = self.get_dim(M)
        return dim[0]==dim[1]
    
    def is_singular(self, M):
        if not self.is_square(M):
            return False
            
        if(self.get_det(M)==0):
            return True
        else:
            return False
    
    def get_transpose(self, M):
        
        if not isinstance(M[0], list):
            M = [M]
            
        dim = self.get_dim(M)
        rows = dim[0]
        cols = dim[1]
        transpose = []
        for k in range(cols):
            col = [M[r][k] for r in range(rows)]
            transpose.append(col)
            
        return transpose
    
    def augment(self, M, b):
        augment = self.copy_matrix(M)
        k = 0
        for i in augment:
            i.append(b[k][0])
            k+=1
        
        return augment
    
    def gaussian_elimination(self, A, b):
        augmented_matrix = self.augment(A, b)
        rows, cols = self.get_dim(augmented_matrix)
        
        for i in range(rows):
            pivotrow = i
            for j in range(i+1, rows):
                if abs(augmented_matrix[j][i]) > abs(augmented_matrix[pivotrow][i]):
                    pivotrow = j
            for k in range(i, cols):
                temp = augmented_matrix[i][k]
                augmented_matrix[i][k] = augmented_matrix[pivotrow][k]
                augmented_matrix[pivotrow][k] = temp
            for j in range(i+1, rows):
                temp = augmented_matrix[j][i]/augmented_matrix[i][i]
                for k in range(i, cols):
                    augmented_matrix[j][k] = float("%.2f" % (augmented_matrix[j][k] - augmented_matrix[i][k] * temp))
        
        return augmented_matrix
    
    def print_gaussian(self, A, b):
        augmented_matrix = self.gaussian_elimination(A, b)
        
        for row in augmented_matrix:
            row = [str(float("%.2f"%i)) for i in row]
            print("\t".join(row))
    
    def get_nature_of_solutions(self, A, b):
        augmented_matrix = self.gaussian_elimination(A, b)
        if(augmented_matrix[-1][-2]==0 and augmented_matrix[-1][-1]==0):
            print("Consistent, singular and infinite solutions")
        elif(augmented_matrix[-1][-2]!=0 and augmented_matrix[-1][-1]==0):
            print("Inconsistent, singular and no solution")
        else:
            print("Consistent, non-singular, unique solution")
    
    def get_upper_triangular(self, A):
        rows, cols = self.get_dim(A)
        C = self.copy_matrix(A)
        
        for i in range(rows):
            pivotrow = i
            for j in range(i+1, rows):
                if abs(C[j][i]) > abs(C[pivotrow][i]):
                    pivotrow = j
            for k in range(i, cols):
                temp = C[i][k]
                C[i][k] = C[pivotrow][k]
                C[pivotrow][k] = temp
            for j in range(i+1, rows):
                temp = C[j][i]/C[i][i]
                for k in range(i, cols):
                    C[j][k] = C[j][k] - C[i][k] * temp
        
        return C
    
    def get_rank(self, A):
        upper = self.get_upper_triangular(A)
        rows, cols = self.get_dim(A)
        
        rank = cols
        for row in upper:
            if(len(set(row))==1 and list(set(row))[0]==0):
                rank -= 1
        return rank
    
    #TODO: Solve system by backward substitution
