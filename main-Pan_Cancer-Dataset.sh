#################### GSE116237 ####################
cat geneList.txt | while read id
do 
echo $id
python ScMutation.py \
    -g $id \
    -s 'ref/GSE116237_Candidate_cancer_expr.h5' \
    -b 'ref/scMutation_ref.h5' \
    -o 'Pan_cancer_Pathyway_detection/GSE116237/Pan_cancer_'$id \
    --lr 1e-4 \
    --batch_size 128 \
    --max_epoch 100 \
    --n_neighbors_pseudo 30 \
    --cl_weight 0.895 \

done

#################### GSE132465 ####################
cat geneList.txt | while read id
do 
echo $id
python ScMutation.py \
    -g $id \
    -s 'ref/GSE132465_Candidate_cancer_expr.h5' \
    -b 'ref/scMutation_ref.h5' \
    -o 'Pan_cancer_Pathyway_detection/GSE132465/Pan_cancer_'$id \
    --lr 1e-4 \
    --batch_size 128 \
    --max_epoch 100 \
    --n_neighbors_pseudo 30 \
    --cl_weight 0.895 \

done

#################### PRJNA591860 ####################
cat geneList.txt | while read id
do 
echo $id
python ScMutation.py \
    -g $id \
    -s 'ref/PRJNA591860_Candidate_cancer_expr.h5' \
    -b 'ref/scMutation_ref.h5' \
    -o 'Pan_cancer_Pathyway_detection/PRJNA591860/Pan_cancer_'$id \
    --lr 1e-4 \
    --batch_size 128 \
    --max_epoch 100 \
    --n_neighbors_pseudo 30 \
    --cl_weight 0.895 \

done

#################### cohort1 ####################
cat geneList.txt | while read id
do 
echo $id
python ScMutation.py \
    -g $id \
    -s 'ref/cohort1_Candidate_cancer_expr.h5' \
    -b 'ref/scMutation_ref.h5' \
    -o 'Pan_cancer_Pathyway_detection/cohort1/Pan_cancer_'$id \
    --lr 1e-4 \
    --batch_size 128 \
    --max_epoch 100 \
    --n_neighbors_pseudo 30 \
    --cl_weight 0.895 \

done

#################### GSE186344 ####################
cat geneList.txt | while read id
do 
echo $id
python ScMutation.py \
    -g $id \
    -s 'ref/GSE186344_Candidate_cancer_expr.h5' \
    -b 'ref/scMutation_ref.h5' \
    -o 'Pan_cancer_Pathyway_detection/GSE186344/Pan_cancer_'$id \
    --lr 1e-4 \
    --batch_size 128 \
    --max_epoch 100 \
    --n_neighbors_pseudo 30 \
    --cl_weight 0.895 \

done


#################### GSE176078 ####################
cat geneList.txt | while read id
do 
echo $id
python ScMutation.py \
    -g $id \
    -s 'ref/GSE176078_Candidate_cancer_expr.h5' \
    -b 'ref/scMutation_ref.h5' \
    -o 'Pan_cancer_Pathyway_detection/GSE176078/Pan_cancer_'$id \
    --lr 1e-4 \
    --batch_size 128 \
    --max_epoch 100 \
    --n_neighbors_pseudo 30 \
    --cl_weight 0.895 \

done
