# Option 1. Creating an environment
Create an environment based on the dependencies.

### Conda environment
Create environment from yaml file:
```
conda env create --file environment.yml
```

Activate:
```
conda activate tbws-ha
```

### Virtual environment
Create environment:
```
python3 -m venv tbws-ha-venv
```

Activate:
```
tbws-ha-venv\Scripts\activate
```

Install requirements: 
```
pip install -r requirements.txt
```

# Option 2. Install packages into existing environment
Or instead of creating an environment, update an existing environment with the dependencies.

### Conda environment
Update environment:
```
conda env update --name env_name --file environment.yml
```

### Virtual environment
Install requirements:
```
pip install -r requirements.txt
```