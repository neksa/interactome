import math
import numpy as np


def mat_to_quat(m):
    """
    Based on http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche52.html
    """
    def sign(x):
        if x >= 0.0: return +1.0
        else: return -1.0

    r11 = m[0, 0]
    r22 = m[1, 1]
    r33 = m[2, 2]
    r12 = m[0, 1]
    r21 = m[1, 0]
    r13 = m[0, 2]
    r31 = m[2, 0]
    r32 = m[2, 1]
    r23 = m[1, 2]

    q0 = (r11 + r22 + r33 + 1.0) / 4.0
    q1 = (r11 - r22 - r33 + 1.0) / 4.0
    q2 = (-r11 + r22 - r33 + 1.0) / 4.0
    q3 = (-r11 - r22 + r33 + 1.0) / 4.0
    if q0 < 0.0: q0 = 0.0
    if q1 < 0.0: q1 = 0.0
    if q2 < 0.0: q2 = 0.0
    if q3 < 0.0: q3 = 0.0
    q0 = math.sqrt(q0)
    q1 = math.sqrt(q1)
    q2 = math.sqrt(q2)
    q3 = math.sqrt(q3)
    if q0 >= q1 and q0 >= q2 and q0 >= q3:
        q0 *= +1.0
        q1 *= sign(r32 - r23)
        q2 *= sign(r13 - r31)
        q3 *= sign(r21 - r12)
    elif q1 >= q0 and q1 >= q2 and q1 >= q3:
        q0 *= sign(r32 - r23)
        q1 *= +1.0
        q2 *= sign(r21 + r12)
        q3 *= sign(r13 + r31)
    elif q2 >= q0 and q2 >= q1 and q2 >= q3:
        q0 *= sign(r13 - r31)
        q1 *= sign(r21 + r12)
        q2 *= +1.0
        q3 *= sign(r32 + r23)
    elif q3 >= q0 and q3 >= q1 and q3 >= q2:
        q0 *= sign(r21 - r12)
        q1 *= sign(r31 + r13)
        q2 *= sign(r32 + r23)
        q3 *= +1.0
    else:
        print "Quaternion: Coding error"
        pass

    q = np.array([q0, q1, q2, q3])
    q /= np.linalg.norm(q)
    return q
