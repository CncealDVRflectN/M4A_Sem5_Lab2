class Matrix(var lineNum: Int, var columnNum: Int) {
    var matrix: Array<Array<Double>> = Array(lineNum, { Array(columnNum, { 0.0 }) })

    constructor(initMtr: Array<Array<Double>>) : this(initMtr.size, initMtr[0].size) {
        for (i in matrix.indices) {
            for (j in matrix[0].indices) {
                matrix[i][j] = initMtr[i][j]
            }
        }
    }

    constructor(initMtr: Matrix) : this(initMtr.matrix)

    fun swapColumns(from: Int, to: Int) {
        var tmp: Double
        for (i in 0 until lineNum) {
            tmp = matrix[i][from]
            matrix[i][from] = matrix[i][to]
            matrix[i][to] = tmp
        }
    }

    fun transpose(): Matrix {
        var result = Matrix(columnNum, lineNum)
        for (i in 0 until columnNum) {
            for (j in 0 until lineNum) {
                result[i][j] = matrix[j][i]
            }
        }
        return result
    }

    fun print() {
        for (i in 0 until lineNum) {
            for (j in 0 until columnNum) {
                print(matrix[i][j].toString() + " ")
            }
            println()
        }
    }

    fun getColumn(j: Int): Array<Double> {
        var result: Array<Double>
        if (j < 0 || j >= columnNum) {
            throw Exception("Incorrect index")
        }
        result = Array(lineNum, { 0.0 })
        for (i in 0 until lineNum) {
            result[i] = matrix[i][j]
        }
        return result
    }

    fun getColumnVector(j: Int): Vector {
        return Vector(getColumn(j))
    }

    fun setColumn(j: Int, column: Array<Double>) {
        if (j < 0 || j >= columnNum) {
            throw Exception("Incorrect index")
        }
        for (i in 0 until lineNum) {
            matrix[i][j] = column[i]
        }
    }

    fun getLine(i: Int): Array<Double> {
        var result: Array<Double>
        if (i < 0 || i >= lineNum) {
            throw Exception("Incorrect index")
        }
        result = Array(columnNum, { 0.0 })
        for (j in 0 until columnNum) {
            result[i] = matrix[i][j]
        }
        return result
    }

    fun setLine(i: Int, line: Array<Double>) {
        if (i < 0 || i >= lineNum) {
            throw Exception("Incorrect index")
        }
        for (j in 0 until columnNum) {
            matrix[i][j] = line[j]
        }
    }

    fun getPart(vararg indices: Int): Matrix {
        var result: Matrix
        if (indices.isEmpty()) {
            throw Exception("Incorrect index")
        }
        result = Matrix(lineNum, indices.size)
        for (i in indices.indices) {
            if (indices[i] < 0 || indices[i] >= columnNum) {
                throw Exception("Incorrect index")
            }
            result.setColumn(i, this.getColumn(indices[i]))
        }
        return result
    }

    fun getPart(indices: Array<Int>): Matrix {
        var result: Matrix
        if (indices.isEmpty()) {
            throw Exception("Incorrect index")
        }
        result = Matrix(lineNum, indices.size)
        for (i in indices.indices) {
            if (indices[i] < 0 || indices[i] >= columnNum) {
                throw Exception("Incorrect index")
            }
            result.setColumn(i, this.getColumn(indices[i]))
        }
        return result
    }

    fun deleteColumns(vararg indices: Int) {
        var newMtr: Array<Array<Double>>
        var counter = 0
        if (indices.size >= columnNum) {
            throw Exception("Incorrect index")
        }
        newMtr = Array(lineNum, { Array(columnNum - indices.size, { 0.0 }) })
        for (j in 0 until columnNum) {
            if (!indices.contains(j)) {
                for (i in 0 until lineNum) {
                    newMtr[i][j - counter] = matrix[i][j]
                }
            } else {
                counter++
            }
        }
        columnNum -= indices.size
        matrix = newMtr
    }

    fun deleteColumns(indices: Array<Int>) {
        var newMtr: Array<Array<Double>>
        var counter = 0
        if (indices.size >= columnNum) {
            throw Exception("Incorrect index")
        }
        newMtr = Array(lineNum, { Array(columnNum - indices.size, { 0.0 }) })
        for (j in 0 until columnNum) {
            if (!indices.contains(j)) {
                for (i in 0 until lineNum) {
                    newMtr[i][j - counter] = matrix[i][j]
                }
            } else {
                counter++
            }
        }
        columnNum -= indices.size
        matrix = newMtr
    }

    operator fun times(b: Vector): Vector {
        var result: Vector
        if (columnNum != b.size) {
            throw Exception("Incorrect size")
        }
        result = Vector(lineNum)
        for (i in 0 until lineNum) {
            result[i] = (0 until columnNum).sumByDouble { matrix[i][it] * b[it] }
        }
        return result
    }

    operator fun get(i: Int): Array<Double> {
        return matrix[i]
    }
}