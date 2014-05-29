(import (rnrs) (rnrs base) (rnrs mutable-pairs))

; cycles: periodic lists ==================================
;===================================================

; create a cyclic list out of its arguments, by pointing
; cdr of the first pair that has () as cdr to the beginning
; of the list. warning: the result is not a proper list.
(define cycle
  (lambda lst
    (letrec ((loop (lambda (first tail)
                  (if (null? (cdr tail))
                    (set-cdr! tail first)
                    (loop first (cdr tail))))))
      (loop lst lst)
      lst)))

; convert a list to a cyclic list
(define list->cycle
  (lambda (lst)
    (apply cycle lst)))

; calculate the length of a cyclic list
(define cycle-length
  (lambda (cycle)
    (let loop ((tail cycle) (n 1))
      (if (eq? (cdr tail) cycle)
        n
        (loop (cdr tail) (+ n 1))))))

; given a cycle and a number, find the number in the cycle
; and return the cycle that has this number as its [car]
(define cycle-find
  (lambda (cycle x)
    (let loop ((tail cycle))
      (if (eq? (car tail) x)
        tail
        (if (eq? (cdr tail) cycle)
          #f
          (loop (cdr tail)))))))

; permutes a single number or a list of numbers by the recipe
; of a single cycle
(define cycle-permute
  (lambda (cycle x)
     (if (list? x)
        (map (lambda (q) (cycle-permute cycle q)) x)
        (let ((result (cycle-find cycle x)))
          (if result
            (cadr result)
            x)))))

; displays a number of elements from a list, used to test the
; working of cycles
(define display-cycle
  (lambda (cycle n)
    (if (eq? n 0)
      '()
      (begin
        (display (car cycle)) (display " ")
        (display-cycle (cdr cycle) (- n 1))))))

; permutations =============================================
;===================================================

; return a function that calculates the permutation given
; by a list representation of the permutation
(define permutation 
  (lambda (permutation-table)
    (letrec ((P (if (vector? permutation-table)
                  permutation-table
                  (list->vector permutation-table)))
             (run (lambda (x)
                    (if (list? x)
                      (map run x)
                      (vector-ref P x)))))
      run)))

; create a permutation-group from a single generator.
; the loop stops when the identity permutation is encountered.
(define permutation-group
  (lambda (permutation identity)
    (let loop ((group (list identity)))
      (let ((next-permutation (permutation (car group))))
        (if (equal? next-permutation identity)
          group
          (loop (cons next-permutation group)))))))
          
; takes to lists representations G, H of permutations and creates
; a new list of permutations of order O(G)xO(H), by applying G to H.
(define permutation-multiply
  (lambda (G H)
    (let loop ((P '()) (g G))
      (if (null? g)
        P
        (loop (append P (map (permutation (car g)) H)) (cdr g))))))

; cyclic groups: permutations from cyclic notation =============
;===================================================

; create a cyclic group by interpreting the elements of the
; input list as cyclic notation
(define cyclic-group
  (lambda (lst)
    (map list->cycle lst)))

; calculate the order of a cyclic-group by taking the
; least common multiple of the lengths of cycles
(define cyclic-group-order
  (lambda (cycle)
    (apply lcm (map cycle-length cycle))))

; permutes a list by the recipe of a cyclic group representation
; can be used to convert a cyclic-group to a permutation
(define cyclic-group-permute
  (lambda (group lst)
    (if (null? group)
      lst
      (cyclic-group-permute (cdr group) (cycle-permute (car group) lst)))))

; set theory functions =========================================
;===================================================

; find all elements that are in only one of the two lists
(define sym-difference
  (lambda (pred a b)
    (let loop ((tail (append a b)) (result '()))
      (if (null? tail)
        result
        (loop (cdr tail)
              (if (and (find (lambda (x) (pred (car tail) x)) a)
                       (find (lambda (x) (pred (car tail) x)) b))
                result
                (cons (car tail) result)))))))

; create all pair-combinations from a list of items
(define combinations
  (lambda (a)
    (let loop ((tail a) (result '()))
      (if (null? (cdr tail))
        result
        (loop (cdr tail) 
              (append result
                      (map (lambda (x) (cons (car tail) x)) (cdr tail))))))))

; removes duplicate entries from a list
(define unique
  (lambda (pred lst)
    (let loop ((tail lst) (result '()))
      (if (null? (cdr tail))
        (reverse (cons (car tail) result))
        (loop (remp (lambda (x) (pred x (car tail))) (cdr tail))
              (cons (car tail) result))))))

; comparator working on unordered pairs. =========================
;===================================================
(define swap-pair
  (lambda (x)
    (cons (cdr x) (car x))))

(define comp-pair
  (lambda (x y)
    (or (equal? x y) (equal? (swap-pair x) y))))

; some utility functions ==========================================
;===================================================

; range creates a list of numbers running from 0 to n-1
(define range
  (lambda (n)
    (let loop ((i n) (lst '()))
      (if (= i 0)
        lst
        (loop (- i 1) (cons (- i 1) lst))))))

; mprint displays all its argument
(define mprint
  (lambda lst
    (do
      ((l lst (cdr l)))
      ((null? l))
      (display (car l)))))

; lfill fills a given string with characters upto a length n
(define lfill
  (lambda (s c n)
    (if (<= n (string-length s))
      s
      (string-append (make-string (- n (string-length s)) c) s))))

; make-encoder takes an associative list and maps keys to values
(define make-encoder
  (lambda (pred codeset)
    (letrec ((run (lambda (x)
                    (if (list? x)
                      (map run x)
                      (let ((result (assp (lambda (q) (pred q x)) codeset)))
                        (if result
                          (cdr result)
                          0))))))
      run)))

; important data and main routines ========================================
;===================================================
; All data related to the cube is using the following labeling:
;    6.----.7
;  4 /___.5| 
;    |   | |
;  2 |___|/ 3
;   0     1 

; there are the three rotation groups associated with the symmetries
; in a oriented cube: an order 2 flip, order 4 rotation around z-axis
; and order 3 rotation around the great diagonal.
(define cycle-F (cyclic-group '((0 7) (3 4) (6 2) (1 5))))
(define cycle-T (cyclic-group '((0 1 3 2) (4 5 7 6))))
(define cycle-R (cyclic-group '((1 2 4) (6 5 3))))

; these are the six fundamental decompositions of the cube
; the cells are denoted by the first point and the clock wise
; orientation of the other three.
(define D-2     '((0 1 2 4) (3 2 1 7) (5 1 4 7) (6 2 7 4) (7 1 4 2)))
(define D-4     '((0 5 1 7) (0 1 3 7) (0 3 2 7) (0 2 6 7) (0 6 4 7) (0 4 5 7)))
(define D-8     '((0 1 3 5) (0 3 7 5) (0 3 6 7) (0 2 6 3) (0 7 6 5) (0 5 6 4)))
(define D-12    '((0 1 7 5) (0 5 7 6) (7 0 6 2) (7 1 0 2) (4 0 6 5) (3 1 7 2)))
(define D-24a   '((0 5 1 7) (0 1 3 7) (0 3 2 7) (0 2 6 7) (0 6 5 7) (0 6 4 5)))
(define D-24b   '((0 1 3 5) (0 3 7 5) (0 7 4 5) (0 7 6 4) (0 7 3 6) (0 3 2 6)))

(define r8 (range 8))

; we can first define the fundamental rotation groups
(define P-2  (permutation-group (permutation 
    (cyclic-group-permute cycle-F r8)) r8))
(define P-3  (permutation-group (permutation 
    (cyclic-group-permute cycle-R r8)) r8))
(define P-4  (permutation-group (permutation 
    (cyclic-group-permute cycle-T r8)) r8))

; and the multiplication groups
(define P-8  (permutation-multiply P-2 P-4))
(define P-12 (permutation-multiply P-4 P-3))
(define P-24 (permutation-multiply P-2 P-12))

; the pressence of an edge is encoded in a bit
; this is not the most efficient method of encoding,
; but it is the easiest.
(define D-edge-encoding
  '(((0 . 3) .    1)
    ((0 . 5) .    2)
    ((0 . 6) .    4)
    ((7 . 1) .    8)
    ((7 . 2) .   16)
    ((7 . 4) .   32)
    ((0 . 7) .   64)
    ((1 . 6) .  128)
    ((2 . 5) .  256)
    ((3 . 4) .  512)))

(define D-latice
  '((0 . 1) (0 . 2) (0 . 4)
    (7 . 5) (7 . 6) (7 . 3)
    (1 . 5) (5 . 4) (4 . 6)
    (6 . 2) (2 . 3) (3 . 1)))

(define D-encode (make-encoder comp-pair D-edge-encoding))

; computes all edges from a set of cells, removing duplicates
; from neighbouring cells
(define cube->edges
  (lambda (cell-list)
    (unique comp-pair 
            (apply append (map combinations cell-list)))))

; given a permutation group and a structure of items,
; we compute this structure under each of the permutations
; in the group
(define permute-dice
  (lambda (P D)
    (map (lambda (p) ((permutation p) D)) P)))

; calculate an id from the binary encoding given above
(define D-id
  (lambda (x)
    (apply + (D-encode (cube->edges x)))))

; find all decompositions of the cube
(define cube-decompositions 
  (append 
    (permute-dice P-2 D-2)
    (permute-dice P-4 D-4)
    (permute-dice P-8 D-8)
    (permute-dice P-12 D-12)
    (permute-dice P-24 D-24a)
    (permute-dice P-24 D-24b)))

(define format-c-list
  (lambda (x)
    (let loop ((tail x) (str ""))
      (if (null? (cdr tail))
        (string-append str (car tail))
        (loop (cdr tail) (string-append 
                           str (car tail) ", "))))))

(define format-cell
  (lambda (x)
    (string-append "Cell(" (format-c-list (map number->string x)) ")")))


; print out
(do
  ((c cube-decompositions (cdr c)) (i 1 (+ i 1)))
  ((null? c))
;  (display (lfill (number->string i) #\  3))
  (display "{ ")
  (display (lfill (number->string (D-id (car c))) #\  4))
  (display ", ")
  (display (length (car c)))
  (display ", { ")
  (display (format-c-list (map format-cell (car c))))
  (display " } },") (newline))

; vim:ts=4:sw=4:et:ci
