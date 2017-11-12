class Vector(var size: Int) {
    var vector: Array<Double> = Array(size, { 0.0 })

    constructor(initVector: Array<Double>) : this(initVector.size) {
        for (i in vector.indices) {
            vector[i] = initVector[i]
        }
    }

    constructor(initVector: Vector) : this(initVector.vector)

    fun swap(from: Int, to: Int) {
        var tmp = vector[from]
        vector[from] = vector[to]
        vector[to] = tmp
    }

    fun print() {
        for (i in vector.indices) {
            print(vector[i].toString() + " ")
        }
        println()
    }

    fun getPartVector(indices: Array<Int>): Vector {
        var result = Vector(indices.size)
        for (i in 0 until indices.size) {
            if (indices[i] < 0 || indices[i] >= size) {
                throw Exception("Incorrect index")
            }
            result[i] = vector[indices[i]]
        }
        return result
    }

    fun getPartVector(indexFrom: Int, indexTo: Int): Vector {
        var result: Vector
        if (indexTo - indexFrom < 0 || indexFrom < 0 || indexFrom >= size || indexTo < 0 || indexFrom >= size) {
            throw Exception("Incorrect indices")
        }
        result = Vector(indexTo - indexFrom + 1)
        for (i in indexFrom..indexTo) {
            result[i - indexFrom] = vector[i]
        }
        return result
    }

    fun min(): Pair<Int, Double> {
        var min = Double.POSITIVE_INFINITY
        var index = -1
        for (i in vector.indices) {
            if (vector[i] < min) {
                min = vector[i]
                index = i
            }
        }
        return Pair(index, min)
    }

    fun max(): Pair<Int, Double> {
        var max = Double.NEGATIVE_INFINITY
        var index = -1
        for (i in vector.indices) {
            if (vector[i] > max) {
                max = vector[i]
                index = i
            }
        }
        return Pair(index, max)
    }

    fun norm(): Double {
        var max = Double.NEGATIVE_INFINITY
        for (i in vector.indices) {
            if (Math.abs(vector[i]) > max) {
                max = Math.abs(vector[i])
            }
        }
        return max
    }

    fun deleteComponents(vararg indices: Int) {
        var newVect: Array<Double>
        var counter = 0
        if (indices.size >= size) {
            throw Exception("Incorrect index")
        }
        newVect = Array(size - indices.size, { 0.0 })
        for (i in 0 until size) {
            if (!indices.contains(i)) {
                newVect[i - counter] = vector[i]
            } else {
                counter++
            }
        }
        size -= indices.size
        vector = newVect
    }

    fun deleteComponents(indices: Array<Int>) {
        var newVect: Array<Double>
        var counter = 0
        if (indices.size >= size) {
            throw Exception("Incorrect index")
        }
        newVect = Array(size - indices.size, { 0.0 })
        for (i in 0 until size) {
            if (!indices.contains(i)) {
                newVect[i - counter] = vector[i]
            } else {
                counter++
            }
        }
        size -= indices.size
        vector = newVect
    }

    operator fun plus(b: Vector): Vector {
        var result: Vector
        if (size != b.size) {
            throw Exception("Incorrect size")
        }
        result = Vector(this)
        for (i in vector.indices) {
            result[i] += b[i]
        }
        return result
    }

    operator fun minus(b: Vector): Vector {
        var result: Vector
        if (size != b.size) {
            throw Exception("Incorrect size")
        }
        result = Vector(this)
        for (i in vector.indices) {
            result[i] -= b[i]
        }
        return result
    }

    operator fun times(b: Double): Vector {
        var result = Vector(this)
        for (i in vector.indices) {
            result[i] *= b
        }
        return result
    }

    operator fun times(b: Vector): Double {
        if (size != b.size) {
            throw Exception("Incorrect size")
        }
        return (0 until size).sumByDouble { vector[it] * b[it] }
    }

    operator fun div(b: Double): Vector {
        var result = Vector(this)
        for (i in vector.indices) {
            result[i] /= b
        }
        return result
    }

    operator fun get(i: Int): Double {
        return vector[i]
    }

    operator fun set(i: Int, b: Double) {
        vector[i] = b
    }
}