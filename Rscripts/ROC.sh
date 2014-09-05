# sort
  sort -r -n b1.scored-label > b1.sorted.scored-label
  sort -r -n b2.scored-label > b2.sorted.scored-label
  sort -r -n b3.scored-label > b3.sorted.scored-label
  sort -r -n b4.scored-label > b4.sorted.scored-label
  sort -r -n b5.scored-label > b5.sorted.scored-label
  sort  -n b6.scored-label > b6.sorted.scored-label
  sort  -n b7.scored-label > b7.sorted.scored-label
  croc-curve -c roc <b1.sorted.scored-label >b1.out
  croc-curve -c roc <b2.sorted.scored-label >b2.out
  croc-curve -c roc <b3.sorted.scored-label >b3.out
  croc-curve -c roc <b4.sorted.scored-label >b4.out
  croc-curve -c roc <b5.sorted.scored-label >b5.out
  croc-curve -c roc <b6.sorted.scored-label >b6.out
  croc-curve -c roc <b7.sorted.scored-label >b7.out

