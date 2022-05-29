# graphslam
Probabilistic Robotics를 참고하여 구현한 GraphSLAM 매트랩 코드입니다.

# Result
파란색 점: motion을 기반으로 그린 input vertex(pose)
빨간색 실선: GraphSLAM 후의 vertex(pose) 궤적
빨간색 점: GraphSLAM 후의 landmark vertex
- 단순한 루프
<p align="center"><img src="/fig/result_1.png"></p>
- 다수의 vertex에 대한 loop
<p align="center"><img src="/fig/result_2.png"></p>