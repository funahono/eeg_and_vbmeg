from itertools import combinations
from collections import Counter
import numpy as np

def predict_ovo_with_tiebreak(modelkun, testkun):
    """
    one versus one での投票数 ＋ 同点時のマージン合計でタイブレークする予測関数

    Parameters
    ----------
    modelkun : sklearn.svm.SVC
        OvOで訓練されたSVMモデル（break_ties=False）
    testkun : array-like, shape (n_samples, n_features)
        テストデータ

    Returns
    -------
    preds : np.ndarray, shape (n_samples,)
        予測ラベル配列
    """
    n_samples = testkun.shape[0]
    classes = modelkun.classes_
    # print(classes)
    n_classes = len(classes)
    class_pairs = list(combinations(classes, 2))  #(classes) Combination 2 の組み合わせをリストとして出してくれる　(ex) classes = 3 → 3C2 :(0,1),(0,2),(1,2)
    # print(class_pairs)

    decision = modelkun.decision_function(testkun)  # shape = (n_samples, n_class_pairs)
    # print(decision)
    preds = []

    for i in range(n_samples):
        # print(f'number of sample:{i}')
        sample_decision = decision[i]
        votes = Counter()
        margin_sum = Counter()

        for idx, (a, b) in enumerate(class_pairs):
            margin = sample_decision[idx]
            winner = a if margin > 0 else b #margin＞0ならa,margin<0ならbに分類、分類されたほうをwinnerとする
            votes[winner] += 1 #votes：分類された回数、winnerは+1ずつ増える
            margin_sum[winner] += abs(margin) #margin_sum：test dataと決定境界の距離の合計値、winnerのみ

        # tei-break：投票数が一番多いクラスを取得
        most_votes = votes.most_common() #(label,votes)のタプルのリストをvotesの多い順に並び替える
        top_vote = most_votes[0][1]#最も出現回数が多いlabelのvotesの数を取得
        top_classes = [cls for cls, v in most_votes if v == top_vote] ##top_votesをvotesとしてもつlabelのリストを取得
        

        # tei-break：votesが同数の場合に、margin_sumが大きい方で評価（より決定関数から距離が大きい方）
        if len(top_classes) == 1: #top_classesに含まれるラベルが1つならそのラベルを取得
            pred = top_classes[0]
        else:#top_classesに含まれるラベルが1つじゃないなら
            pred = max(top_classes, key=lambda cls: margin_sum[cls]) #各ラベルのmargin_sumを取得して比較、大きいクラスを選択する
        
        # print(int(pred))

        preds.append(pred)
    
    # print(decision)

    return np.array(preds), decision
