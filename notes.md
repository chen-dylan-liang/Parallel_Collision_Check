- P1

  Hi, I'm Dylan. I'm doing this final project with Qiang. Parallel 3D Collision Checking (Written in Rust) 

- P2

  To begin with, I wanna introduce the research significance of collision checking and why it makes sense to apply parallel computing to accelerate the collision checking algorithm.

  

  Here are pictures of the robot dog produced by Boston Dynamics. The key algorithm to guarantee the safety of this robot to make sure it doesn't hit the ladders, is collision checking. In these applications, computational speed is of the main concerns. We wanna make it real-time, to make the robot locomotions safe but agile. So here comes parallel computing. 

- P3

  Collision checking algorithms involves two phases, a broad one and a narrow one. In broad phase, a tree structure called BVH is built and traversed. 3D BVH is similar to a kd-tree with k=3, but it splits the space not in terms of actual positions, but in terms of groups of 3D primitives. By traversing BVH, primitives far away, without any possibilities to collide, are ruled out, leaving potential colliding pairs to the narrow phase. 

- P4

  In narrow phase, a geometric algorithm called GJK will be run across all pairs to decided whether they really collide.

- P5

  The process of doing this algorithm in parallel is much similar to homework 5. First, build a tree structure, BVH, in parallel, next, traverse the tree in parallel to extract pairs to narrow phase. Then narrow phase can be embarrassingly parallelized, just like how we do queries on a kd-tree in parallel. However, the details in this algorithm are actually far more delicate than the ones in kd-tree. For example, in the traversing part, because we need to traverse pairs, so we're traversing one tree and its own copy in parallel instead of just traversing the tree itself.

- P6

  Finally, I wanna justify our usage of rust rather than C++. Rust is as fast as C++, but compared to C++, rust has friendlier dependency management mechanisms, so the parallel computing library it offers has the ability to be migrated to different machines and operating systems with zero pain. It makes the workflow easier.

  











