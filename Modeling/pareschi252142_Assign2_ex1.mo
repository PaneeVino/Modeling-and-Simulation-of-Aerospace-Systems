within ;
model pareschi252142_Assign2_ex1
model Thermal_model

model Var_eps
    Modelica.Blocks.Interfaces.RealInput theta
  annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=270,
        origin={0,100}),iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=270,
        origin={0,100})));
     Modelica.Blocks.Interfaces.RealOutput eps
  annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=270,
        origin={0,-100}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=270,
        origin={0,-100})));
  parameter Real eps_min(  unit="-");
  parameter Real eps_max(  unit="-");
equation
  eps = eps_min + (eps_max-eps_min)/(0.4*Modelica.Constants.pi)*(theta+0.4*Modelica.Constants.pi);
  connect(eps, eps) annotation (Line(points={{0,-100},{0,-100}},
                                                             color={0,0,127}));
  annotation (Diagram(graphics={Text(
          extent={{36,30},{-38,-8}},
          textColor={28,108,200},
          textString="Eps_Var")}));
end Var_eps;

model Radiator_variable
  extends Modelica.Thermal.HeatTransfer.Interfaces.Element1D;
  Modelica.Blocks.Interfaces.RealInput eps
    annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=270,
        origin={2,98}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=270,
        origin={2,98})));
        parameter Real A(  unit="m2");

equation
  Q_flow = eps*A*Modelica.Constants.sigma*(port_a.T^4 - port_b.T^4);
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
            100,100}}), graphics={
        Rectangle(
          extent={{50,80},{90,-80}},
          fillColor={192,192,192},
          fillPattern=FillPattern.Backward),
        Rectangle(
          extent={{-90,80},{-50,-80}},
          fillColor={192,192,192},
          fillPattern=FillPattern.Backward),
        Line(points={{-36,10},{36,10}}, color={191,0,0}),
        Line(points={{-36,10},{-26,16}}, color={191,0,0}),
        Line(points={{-36,10},{-26,4}}, color={191,0,0}),
        Line(points={{-36,-10},{36,-10}}, color={191,0,0}),
        Line(points={{26,-16},{36,-10}}, color={191,0,0}),
        Line(points={{26,-4},{36,-10}}, color={191,0,0}),
        Line(points={{-36,-30},{36,-30}}, color={191,0,0}),
        Line(points={{-36,-30},{-26,-24}}, color={191,0,0}),
        Line(points={{-36,-30},{-26,-36}}, color={191,0,0}),
        Line(points={{-36,30},{36,30}}, color={191,0,0}),
        Line(points={{26,24},{36,30}}, color={191,0,0}),
        Line(points={{26,36},{36,30}}, color={191,0,0}),
        Text(
          extent={{-150,125},{150,85}},
          textString="%name",
          textColor={0,0,255}),
        Text(
          extent={{-150,-90},{150,-120}},
          textString="Gr=%Gr"),
        Rectangle(
          extent={{-50,80},{-44,-80}},
          lineColor={191,0,0},
          fillColor={191,0,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{45,80},{50,-80}},
          lineColor={191,0,0},
          fillColor={191,0,0},
          fillPattern=FillPattern.Solid)}),
    Documentation(info="<html>
<p>
This is a model describing the thermal radiation, i.e., electromagnetic
radiation emitted between two bodies as a result of their temperatures.
The following constitutive equation is used:
</p>
<blockquote><pre>
Q_flow = Gr*sigma*(port_a.T^4 - port_b.T^4);
</pre></blockquote>
<p>
where Gr is the radiation conductance and sigma is the Stefan-Boltzmann
constant (= Modelica.Constants.sigma). Gr may be determined by
measurements and is assumed to be constant over the range of operations.
</p>
<p>
For simple cases, Gr may be analytically computed. The analytical
equations use epsilon, the emission value of a body which is in the
range 0..1. Epsilon=1, if the body absorbs all radiation (= black body).
Epsilon=0, if the body reflects all radiation and does not absorb any.
</p>
<blockquote><pre>
Typical values for epsilon:
aluminium, polished    0.04
copper, polished       0.04
gold, polished         0.02
paper                  0.09
rubber                 0.95
silver, polished       0.02
wood                   0.85..0.9
</pre></blockquote>
<p><strong>Analytical Equations for Gr</strong></p>
<p>
<strong>Small convex object in large enclosure</strong>
(e.g., a hot machine in a room):
</p>
<blockquote><pre>
Gr = e*A
where
   e: Emission value of object (0..1)
   A: Surface area of object where radiation
      heat transfer takes place
</pre></blockquote>
<p><strong>Two parallel plates</strong>:</p>
<blockquote><pre>
Gr = A/(1/e1 + 1/e2 - 1)
where
   e1: Emission value of plate1 (0..1)
   e2: Emission value of plate2 (0..1)
   A : Area of plate1 (= area of plate2)
</pre></blockquote>
<p><strong>Two long cylinders in each other</strong>, where radiation takes
place from the inner to the outer cylinder):
</p>
<blockquote><pre>
Gr = 2*pi*r1*L/(1/e1 + (1/e2 - 1)*(r1/r2))
where
   pi: = Modelica.Constants.pi
   r1: Radius of inner cylinder
   r2: Radius of outer cylinder
   L : Length of the two cylinders
   e1: Emission value of inner cylinder (0..1)
   e2: Emission value of outer cylinder (0..1)
</pre></blockquote>
</html>"),
    uses(Modelica(version="4.0.0")));
end Radiator_variable;

  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature DeepSpace(T(
        displayUnit="K") = 3) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={0,-68})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor MainBody(C=1.5*10^5, T(
      fixed=true,
      start=298.15,
      displayUnit="K"))
    annotation (Placement(transformation(extent={{-14,70},{12,96}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Panel2(C=1187.5, T(
      fixed=true,
      start=298.15,
      displayUnit="K"))
    annotation (Placement(transformation(extent={{68,74},{94,100}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Panel1(C=1187.5, T(
      fixed=true,
      start=298.15,
      displayUnit="K"))
    annotation (Placement(transformation(extent={{-98,74},{-72,100}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor G12(G=10)
    annotation (Placement(transformation(extent={{-40,56},{-20,76}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor G13(G=10)
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={30,66})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Radiator1(C=30, T(
      fixed=true,
      start=298.15,
      displayUnit="K")) annotation (Placement(transformation(
        extent={{-12,-12},{12,12}},
        rotation=90,
        origin={-112,0})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Radiator2(C=30, T(
      fixed=true,
      start=298.15,
      displayUnit="K")) annotation (Placement(transformation(
        extent={{-12,-12},{12,12}},
        rotation=270,
        origin={112,0})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor G14(G=10)
    annotation (Placement(transformation(extent={{-40,6},{-20,26}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor G15(G=10)
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={32,16})));
  Modelica.Thermal.HeatTransfer.Components.BodyRadiation eps1A1(Gr=0.45*3.5)
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={0,-8})));
  Modelica.Thermal.HeatTransfer.Components.BodyRadiation eps3A3(Gr=0.75*0.475)
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={68,42})));
  Modelica.Thermal.HeatTransfer.Components.BodyRadiation eps2A2(Gr=0.75*0.475)
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-64,46})));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow alpha2A2Psun(Q_flow=0.78*
        0.475*1350) annotation (Placement(transformation(
        extent={{-14,-14},{14,14}},
        rotation=0,
        origin={-116,66})));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow alpha3A3Psun(Q_flow=0.78*
        0.475*1350) annotation (Placement(transformation(
        extent={{-14,-14},{14,14}},
        rotation=180,
        origin={116,66})));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow alpha1A1topPsun(Q_flow=
        0.25*0.6*1350)
    annotation (Placement(transformation(extent={{-40,28},{-14,54}})));
  Var_eps var_eps(eps_min = 0.01, eps_max = 0.98)  annotation (Placement(transformation(
        extent={{-12,-12},{12,12}},
        rotation=180,
        origin={0,-102})));
  Radiator_variable epsA5(A=0.25) annotation (Placement(transformation(
        extent={{-12,-12},{12,12}},
        rotation=180,
        origin={36,-28})));
  Radiator_variable epsA4(A=0.25) annotation (Placement(transformation(
        extent={{-12,-12},{12,12}},
        rotation=0,
        origin={-32,-28})));
  Modelica.Blocks.Interfaces.RealInput theta annotation (Placement(
        transformation(
        extent={{-7,-7},{7,7}},
        rotation=90,
        origin={-1,-117}), iconTransformation(
        extent={{-7,-7},{7,7}},
        rotation=90,
        origin={-1,-117})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-34,104})));
  Modelica.Blocks.Interfaces.RealOutput y annotation (Placement(transformation(
        extent={{-8,-8},{8,8}},
        rotation=90,
        origin={-26,130}), iconTransformation(
        extent={{-8,-8},{8,8}},
        rotation=90,
        origin={-26,130})));
equation
  connect(G12.port_b, MainBody.port)
    annotation (Line(points={{-20,66},{-1,66},{-1,70}}, color={191,0,0}));
  connect(G12.port_a, Panel1.port)
    annotation (Line(points={{-40,66},{-85,66},{-85,74}}, color={191,0,0}));
  connect(G13.port_a, Panel2.port)
    annotation (Line(points={{40,66},{81,66},{81,74}}, color={191,0,0}));
  connect(G13.port_b, MainBody.port)
    annotation (Line(points={{20,66},{-1,66},{-1,70}}, color={191,0,0}));
  connect(Radiator2.port, G15.port_a) annotation (Line(points={{100,0},{58,0},{
          58,16},{42,16}}, color={191,0,0}));
  connect(G15.port_b, MainBody.port)
    annotation (Line(points={{22,16},{-1,16},{-1,70}}, color={191,0,0}));
  connect(G14.port_b, MainBody.port)
    annotation (Line(points={{-20,16},{-1,16},{-1,70}}, color={191,0,0}));
  connect(Radiator1.port, G14.port_a) annotation (Line(points={{-100,0},{-54,0},
          {-54,16},{-40,16}}, color={191,0,0}));
  connect(eps1A1.port_a, MainBody.port)
    annotation (Line(points={{0,2},{-1,4},{-1,70}}, color={191,0,0}));
  connect(eps1A1.port_b, DeepSpace.port)
    annotation (Line(points={{0,-18},{0,-58}}, color={191,0,0}));
  connect(eps3A3.port_a, Panel2.port) annotation (Line(points={{68,52},{68,66},
          {81,66},{81,74}}, color={191,0,0}));
  connect(eps3A3.port_b, DeepSpace.port) annotation (Line(points={{68,32},{68,
          -52},{0,-52},{0,-58}}, color={191,0,0}));
  connect(eps2A2.port_a, Panel1.port) annotation (Line(points={{-64,56},{-64,66},
          {-85,66},{-85,74}}, color={191,0,0}));
  connect(eps2A2.port_b, DeepSpace.port) annotation (Line(points={{-64,36},{-64,
          -52},{0,-52},{0,-58}}, color={191,0,0}));
  connect(alpha2A2Psun.port, Panel1.port)
    annotation (Line(points={{-102,66},{-85,66},{-85,74}}, color={191,0,0}));
  connect(alpha3A3Psun.port, Panel2.port)
    annotation (Line(points={{102,66},{81,66},{81,74}}, color={191,0,0}));
  connect(alpha1A1topPsun.port, MainBody.port)
    annotation (Line(points={{-14,41},{-1,41},{-1,70}}, color={191,0,0}));
  connect(epsA5.port_a, Radiator2.port) annotation (Line(points={{48,-28},{58,
          -28},{58,0},{100,0}}, color={191,0,0}));
  connect(epsA5.port_b, DeepSpace.port)
    annotation (Line(points={{24,-28},{0,-28},{0,-58}}, color={191,0,0}));
  connect(epsA5.eps, var_eps.eps) annotation (Line(points={{35.76,-39.76},{
          35.76,-82},{0,-82},{0,-90}}, color={0,0,127}));
  connect(epsA4.port_b, DeepSpace.port)
    annotation (Line(points={{-20,-28},{0,-28},{0,-58}}, color={191,0,0}));
  connect(epsA4.port_a, Radiator1.port) annotation (Line(points={{-44,-28},{-54,
          -28},{-54,0},{-100,0}}, color={191,0,0}));
  connect(epsA4.eps, var_eps.eps) annotation (Line(points={{-31.76,-16.24},{
          -31.76,-8},{-72,-8},{-72,-82},{0,-82},{0,-90}}, color={0,0,127}));
  connect(theta, var_eps.theta) annotation (Line(points={{-1,-117},{-4,-117},{
          -4,-114},{0,-114}}, color={0,0,127}));
  connect(temperatureSensor.port, MainBody.port) annotation (Line(points={{-34,
          94},{-34,86},{-14,86},{-14,70},{-1,70}}, color={191,0,0}));
  connect(temperatureSensor.T, y)
    annotation (Line(points={{-34,115},{-26,115},{-26,130}}, color={0,0,127}));
  annotation (
    uses(Modelica(version="4.0.0")),
    Diagram(coordinateSystem(extent={{-140,-120},{140,140}})),
    Icon(coordinateSystem(extent={{-140,-120},{140,140}}), graphics={
        Rectangle(
          extent={{-140,140},{140,-140}},
          lineColor={28,108,200},
          fillColor={244,125,35},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-114,92},{114,-80}},
          textColor={28,108,200},
          textString="Thermal Model"),
        Text(
          extent={{10,-66},{86,-132}},
          textColor={28,108,200},
          textString="theta"),
        Text(
          extent={{-14,140},{40,106}},
          textColor={28,108,200},
          textString="T1")}));
end Thermal_model;

model Control_Logic
  Modelica.Blocks.Interfaces.RealInput theta
    annotation (Placement(transformation(extent={{-140,10},{-100,50}}),
        iconTransformation(extent={{-140,10},{-100,50}})));
  Modelica.Blocks.Interfaces.RealInput T1
    annotation (Placement(transformation(extent={{-140,-88},{-100,-48}}),
        iconTransformation(extent={{-140,-88},{-100,-48}})));
  Modelica.Blocks.Interfaces.RealOutput V annotation (Placement(transformation(
          extent={{100,-30},{140,10}}), iconTransformation(extent={{100,-30},{
            140,10}})));
  Modelica.Blocks.Logical.LogicalSwitch logicalSwitch
    annotation (Placement(transformation(extent={{0,-20},{20,0}})));
  Modelica.Blocks.Logical.GreaterEqualThreshold greaterEqualThreshold(threshold
      =-0.4*Modelica.Constants.pi)
    annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
  Modelica.Blocks.Logical.And and1
    annotation (Placement(transformation(extent={{-40,0},{-20,20}})));
  Modelica.Blocks.Logical.LessEqualThreshold lessEqualThreshold
    annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));
  Modelica.Blocks.Logical.LessEqual lessEqual
    annotation (Placement(transformation(extent={{-40,-40},{-20,-20}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T(
        displayUnit="K") = 294.15)
    annotation (Placement(transformation(extent={{-100,-60},{-80,-40}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
    annotation (Placement(transformation(extent={{-68,-60},{-48,-40}})));
  Modelica.Blocks.Math.Add add(k2=-1)
    annotation (Placement(transformation(extent={{-20,-80},{0,-60}})));
  Modelica.Blocks.Math.Gain gain(k=6*10^(-5))
    annotation (Placement(transformation(extent={{20,-80},{40,-60}})));
  Modelica.Blocks.Logical.Switch switch1
    annotation (Placement(transformation(extent={{60,-20},{80,0}})));
  Modelica.Blocks.Sources.Constant const(k=0)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=180,
        origin={70,-50})));
equation
  connect(greaterEqualThreshold.u, theta) annotation (Line(points={{-82,30},{
          -120,30}},                   color={0,0,127}));
  connect(lessEqualThreshold.u, theta) annotation (Line(points={{-82,-10},{-92,
          -10},{-92,30},{-120,30}},color={0,0,127}));
  connect(lessEqualThreshold.y, and1.u2) annotation (Line(points={{-59,-10},{
          -50,-10},{-50,2},{-42,2}},
                                  color={255,0,255}));
  connect(and1.u1, greaterEqualThreshold.y) annotation (Line(points={{-42,10},{
          -54,10},{-54,30},{-59,30}}, color={255,0,255}));
  connect(and1.y, logicalSwitch.u2) annotation (Line(points={{-19,10},{-16,10},
          {-16,-2},{-8,-2},{-8,-10},{-2,-10}},
                                             color={255,0,255}));
  connect(logicalSwitch.u1, and1.y) annotation (Line(points={{-2,-2},{-16,-2},{
          -16,10},{-19,10}}, color={255,0,255}));
  connect(lessEqual.y, logicalSwitch.u3) annotation (Line(points={{-19,-30},{-8,
          -30},{-8,-18},{-2,-18}},
                               color={255,0,255}));
  connect(temperatureSensor.port, fixedTemperature.port)
    annotation (Line(points={{-68,-50},{-80,-50}}, color={191,0,0}));
  connect(temperatureSensor.T, lessEqual.u2) annotation (Line(points={{-47,-50},
          {-44,-50},{-44,-38},{-42,-38}}, color={0,0,127}));
  connect(add.u2, temperatureSensor.T) annotation (Line(points={{-22,-76},{-42,
          -76},{-42,-50},{-47,-50}}, color={0,0,127}));
  connect(add.u1, T1) annotation (Line(points={{-22,-64},{-70,-64},{-70,-68},{
          -120,-68}}, color={0,0,127}));
  connect(gain.u, add.y)
    annotation (Line(points={{18,-70},{1,-70}},  color={0,0,127}));
  connect(logicalSwitch.y, switch1.u2)
    annotation (Line(points={{21,-10},{58,-10}},
                                               color={255,0,255}));
  connect(gain.y, switch1.u1) annotation (Line(points={{41,-70},{46,-70},{46,-2},
          {58,-2}}, color={0,0,127}));
  connect(const.y, switch1.u3) annotation (Line(points={{59,-50},{50,-50},{50,
          -18},{58,-18}},
                   color={0,0,127}));
  connect(switch1.y, V)
    annotation (Line(points={{81,-10},{120,-10}},color={0,0,127}));
  connect(lessEqual.u1, T1) annotation (Line(points={{-42,-30},{-70,-30},{-70,
          -68},{-120,-68}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,60}}), graphics={
        Rectangle(
          extent={{-100,60},{100,-100}},
          lineColor={28,108,200},
          fillColor={85,255,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-108,56},{-16,26}},
          textColor={28,108,200},
          textString="theta"),
        Text(
          extent={{-106,-62},{-44,-100}},
          textColor={28,108,200},
          textString="T1"),
        Text(
          extent={{68,48},{98,6}},
          textColor={28,108,200},
          textString="V"),
        Text(
          extent={{-96,34},{88,-62}},
          textColor={28,108,200},
          textString="Control Logic")}),                         Diagram(
        coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,60}})),
    uses(Modelica(version="4.0.0")));
end Control_Logic;

model Electro_Mechanical_Model
  Modelica.Electrical.Analog.Basic.Resistor resistor(R=0.1)
    annotation (Placement(transformation(extent={{-24,34},{-4,54}})));
  Modelica.Electrical.Analog.Basic.Ground ground
    annotation (Placement(transformation(extent={{-10,-40},{10,-20}})));
  Modelica.Electrical.Analog.Basic.Inductor inductor(L=0.001)
    annotation (Placement(transformation(extent={{8,34},{28,54}})));
  Modelica.Electrical.Analog.Sources.SignalVoltage signalVoltage annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-44,16})));
  Modelica.Electrical.Analog.Basic.RotationalEMF emf(k=0.3)
    annotation (Placement(transformation(extent={{32,12},{52,32}})));
  Modelica.Mechanics.Rotational.Sensors.AngleSensor angleSensor
    annotation (Placement(transformation(extent={{96,12},{116,32}})));
  Modelica.Mechanics.Rotational.Components.Inertia inertia(
    J=1/60,
    phi(start=-0.4*Modelica.Constants.pi),
    w(start=0)) annotation (Placement(transformation(extent={{64,12},{84,32}})));
  Modelica.Blocks.Interfaces.RealOutput theta
    annotation (Placement(transformation(extent={{120,12},{140,32}})));
  Modelica.Blocks.Interfaces.RealInput V annotation (Placement(transformation(
        extent={{-12,-12},{12,12}},
        rotation=90,
        origin={-20,-30}), iconTransformation(
        extent={{-12,-12},{12,12}},
        rotation=90,
        origin={-20,-30})));
equation
  connect(resistor.n, inductor.p)
    annotation (Line(points={{-4,44},{8,44}}, color={0,0,255}));
  connect(emf.p, inductor.n)
    annotation (Line(points={{42,32},{42,44},{28,44}}, color={0,0,255}));
  connect(ground.p, signalVoltage.n) annotation (Line(points={{0,-20},{0,-14},{
          -44,-14},{-44,6}},color={0,0,255}));
  connect(emf.n, ground.p) annotation (Line(points={{42,12},{42,-14},{0,-14},{0,
          -20}}, color={0,0,255}));
  connect(emf.flange, inertia.flange_a)
    annotation (Line(points={{52,22},{64,22}}, color={0,0,0}));
  connect(inertia.flange_b, angleSensor.flange)
    annotation (Line(points={{84,22},{96,22}}, color={0,0,0}));
  connect(angleSensor.phi, theta)
    annotation (Line(points={{117,22},{130,22}}, color={0,0,127}));
  connect(signalVoltage.v, V)
    annotation (Line(points={{-32,16},{-20,16},{-20,-30}}, color={0,0,127}));
  connect(resistor.p, signalVoltage.p)
    annotation (Line(points={{-24,44},{-44,44},{-44,26}}, color={0,0,255}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-60,-40},{120,60}}),
        graphics={
        Rectangle(
          extent={{-60,60},{120,-40}},
          lineColor={28,108,200},
          fillColor={128,255,0},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{80,60},{118,28}},
          textColor={28,108,200},
          textString="theta"),
        Text(
          extent={{-54,-12},{-34,-36}},
          textColor={28,108,200},
          textString="V"),
        Text(
          extent={{-52,40},{112,-12}},
          textColor={28,108,200},
          textString="Electro-mechanical model")}),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-60,-40},{120,
            60}})),
    uses(Modelica(version="4.0.0")));
end Electro_Mechanical_Model;

  Thermal_model thermal_model annotation (Placement(transformation(
        extent={{-45,-42},{45,42}},
        rotation=90,
        origin={-62,-39})));
  Control_Logic control_Logic
    annotation (Placement(transformation(extent={{-90,34},{-12,92}})));
  Electro_Mechanical_Model electro_Mechanical_Model annotation (Placement(
        transformation(
        extent={{-51,-37},{51,37}},
        rotation=270,
        origin={73,39})));
equation
  connect(thermal_model.y, control_Logic.T1) annotation (Line(points={{-100.769,
          -47.3571},{-112,-47.3571},{-112,45.6},{-97.8,45.6}},          color={
          0,0,127}));
  connect(electro_Mechanical_Model.theta, thermal_model.theta) annotation (Line(
        points={{81.88,-17.6667},{81.88,-39.3214},{-20.9692,-39.3214}}, color={
          0,0,127}));
  connect(control_Logic.theta, electro_Mechanical_Model.theta) annotation (Line(
        points={{-97.8,81.125},{-128,81.125},{-128,-90},{81.88,-90},{81.88,
          -17.6667}}, color={0,0,127}));
  connect(control_Logic.V, electro_Mechanical_Model.V) annotation (Line(points={{-4.2,
          66.625},{36,66.625},{36,67.3333},{43.4,67.3333}},        color={0,0,
          127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-140,
            -100},{120,100}})),
  Diagram(coordinateSystem(preserveAspectRatio = false, extent={{-140,-100},{
            120,100}})),
  version = "",
  uses);
end pareschi252142_Assign2_ex1;
