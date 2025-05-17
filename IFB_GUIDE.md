# Guide d'utilisation du cluster IFB pour l'analyse RNA-seq poulet

## 1. Connexion au cluster et configuration

1. Se connecter au cluster :
```bash
ssh hnguyen97@core.cluster.france-bioinformatique.fr
```

2. Démarrer une session tmux (protection contre les déconnexions) :
```bash
# Créer une nouvelle session tmux
tmux new -s rnaseq

# Si vous êtes déconnecté, vous pouvez revenir à votre session avec :
tmux attach -t rnaseq
```

## 2. Transfert des fichiers depuis votre machine locale

1. Transférer les scripts et données :
```bash
cd /home/nguyeho3/Documents/RNA_seq_test_for_github_repo

# Créer une archive des données nécessaires
tar czf chicken_data.tar.gz \
    chicken_atlas/scripts/download_reference.sh \
    chicken_atlas/genome/* \
    chicken_atlas/SRR13302173/* \
    rna_seq_ifb.sh \
    download_reference_ifb.sh

# Transférer l'archive vers le cluster
scp chicken_data.tar.gz hnguyen97@core.cluster.france-bioinformatique.fr:~/ondemand/data/sys/dashboard/batch_connect/sys/jupyter/core/
```

## 3. Préparation sur le cluster

Dans votre session tmux, exécutez :

1. Créer et configurer le répertoire de travail :
```bash
# Aller dans le répertoire de travail
cd ~/ondemand/data/sys/dashboard/batch_connect/sys/jupyter/core
mkdir -p chicken_rnaseq/{data,genome,scripts,status}
cd chicken_rnaseq

# Extraire les fichiers
tar xzf ../chicken_data.tar.gz

# Organiser les fichiers
mv chicken_atlas/genome/* genome/
mv chicken_atlas/SRR13302173/* data/
mv chicken_atlas/scripts/download_reference.sh scripts/
mv rna_seq_ifb.sh download_reference_ifb.sh ./

# Nettoyer
rm -rf chicken_atlas
rm ../chicken_data.tar.gz

# Rendre les scripts exécutables
chmod +x *.sh scripts/*.sh
```

2. Vérifier les fichiers :
```bash
# Vérifier la structure
ls -R

# Vérifier l'espace disponible
df -h .
```

3. Mettre à jour les chemins dans les scripts :
```bash
# Mettre à jour le chemin de base dans les scripts
WORK_DIR="$PWD"
sed -i "s|WORK_DIR=.*|WORK_DIR=\"${WORK_DIR}\"|" download_reference_ifb.sh
sed -i "s|WORK_DIR=.*|WORK_DIR=\"${WORK_DIR}\"|" rna_seq_ifb.sh
```

4. Télécharger et indexer le génome de référence :
```bash
# Lancer le job avec plus de ressources
sbatch download_reference_ifb.sh

# Surveiller le progrès
watch -n 10 'squeue -u hnguyen97'  # Appuyez sur Ctrl-C pour quitter watch
tail -f download_ref_*.out         # Appuyez sur Ctrl-C pour quitter tail

# En cas de déconnexion pendant l'exécution :
# 1. Reconnectez-vous au cluster
# 2. Rattachez-vous à votre session tmux : tmux attach -t rnaseq
# 3. Le script continuera grâce aux points de reprise
```

5. Une fois le génome prêt, lancer l'analyse :
```bash
# Vérifier que la préparation du génome est terminée
if [ -f "genome/reference_preparation_complete" ]; then
    sbatch rna_seq_ifb.sh
else
    echo "La préparation du génome n'est pas terminée"
fi
```

## 4. Surveillance de l'analyse

1. Dans votre session tmux, vérifier l'état des jobs :
```bash
# Voir tous vos jobs
squeue -u hnguyen97

# Voir les détails d'un job en cours
scontrol show job <JOB_ID>
```

2. Surveiller les logs :
```bash
# Dans un panneau tmux séparé (Ctrl-b puis ")
tail -f rnaseq_*.out
tail -f rnaseq_*.err
```

3. Vérifier l'utilisation des ressources :
```bash
sstat -j <JOB_ID>
```

## 5. Récupération des résultats

Une fois l'analyse terminée, depuis votre machine locale :
```bash
cd /home/nguyeho3/Documents/RNA_seq_test_for_github_repo
scp -r hnguyen97@core.cluster.france-bioinformatique.fr:~/ondemand/data/sys/dashboard/batch_connect/sys/jupyter/core/chicken_rnaseq/results/ ./
```

## Notes importantes

- Le script utilise :
  - La dernière version du génome poulet (GRCg7b) d'Ensembl release 110
  - Les annotations GTF d'Ensembl release 110
  - Points de reprise pour la préparation du génome
  - Protection contre les déconnexions avec tmux
  - 8 CPUs et 32GB de RAM pour l'indexation
- Structure des répertoires sur le cluster :
  - `genome/` : Fichiers du génome et index
  - `data/` : Données RNA-seq
  - `scripts/` : Scripts utilitaires
  - `results/` : Résultats de l'analyse
  - `status/` : Points de reprise
- Commandes tmux utiles :
  - `Ctrl-b d` : Détacher la session
  - `Ctrl-b "` : Diviser horizontalement
  - `Ctrl-b %` : Diviser verticalement
  - `Ctrl-b flèches` : Naviguer entre les panneaux
- En cas de déconnexion :
  1. Reconnectez-vous au cluster
  2. `tmux attach -t rnaseq`
  3. Les jobs SLURM continuent même si vous êtes déconnecté

## Ressources utiles

- Documentation IFB : https://ifb-elixirfr.gitlab.io/cluster/doc/
- État du cluster : https://www.france-bioinformatique.fr/cluster-status/
- Support : https://support.cluster.france-bioinformatique.fr/
- Documentation tmux : https://tmuxcheatsheet.com/
- Génome de référence : https://www.ensembl.org/Gallus_gallus/Info/Index 

notes : 
scp -C /home/nguyeho3/Documents/RNA_seq_test_for_github_repo/chicken_atlas/data/raw/* hnguyen97@core.cluster.france-bioinformatique.fr:~/ondemand/data/sys/dashboard/batch_connect/sys/jupyter/core/chicken_rnaseq/data/raw/