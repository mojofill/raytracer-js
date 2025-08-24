const canvas = document.getElementById('canvas');
const ctx = canvas.getContext('2d');

const WIDTH = 800;
const HEIGHT = 800;

canvas.width = WIDTH;
canvas.height = HEIGHT;

ctx.clearRect(0, 0, canvas.width, canvas.height);
ctx.fillStyle = 'black';
ctx.fillRect(0, 0, canvas.width, canvas.height);

const fov = 110*Math.PI/180;

// opengl view plane: x right, y up, right hand rule z out the page. forward = into page = -3
const camera = {
    position: vec3(0, 0, 10),
    forward: vec3(0, 0, -1),
    right: vec3(1, 0, 0),
    up: vec3(0, 1, 0)
};

const MAX_DEPTH = 100;
const backgroundColor = vec3(0.5, 0.5, 0.5);

const ambientLight = vec3(1.0, 1.0, 1.0);

const lights = [];
lights.push(
    {position: vec3(0, 5, 0), color: vec3(1.0, 1.0, 1.0), radius: 0.5, mat: createMaterial(vec3(0,0,0), 0, vec3(1,1,1))},
    {position: vec3(0, -5, 0), color: vec3(1.0, 1.0, 1.0), radius: 0.5, mat: createMaterial(vec3(0,0,0), 0, vec3(1,1,1))},
    // {position: vec3(-5, 0, 0), color: vec3(1.0, 1.0, 1.0), radius: 0.5, mat: createMaterial(vec3(0,0,0), 0, vec3(1,1,1))},
    // {position: vec3(5, 0, 0), color: vec3(1.0, 1.0, 1.0), radius: 0.5, mat: createMaterial(vec3(0,0,0), 0, vec3(1,1,1))}
);

// let width = 8;
// let num_across = 1;
// let radius = width/num_across/2;

// for (let y = -width/2; y <= width/2; y += 2*radius) {
//     for (let x = -width/2; x <= width/2; x += 2*radius) {
//         objs.push({
//             position: vec3(x,y,0),
//             radius: radius,
//             mat: createMaterial(vec3(0.5,0.5,1),1)
//         })
//     }
// }

let t = 0.9;
let a = 0.7;

const H = 10;

const objs = [];
objs.push(
    {
        position: vec3(0, 0, 0),
        radius: 4.25,
        mat: createMaterial(vec3(0.5,0.5,0.5), 0.9)
    }
)

function generateUVSphere(radius, stacks, slices) {
    let verts = []; // flat array of triangles: [x,y,z, x,y,z, ...]

    for (let i = 0; i < stacks; i++) {
        let theta1 = Math.PI * i / stacks;
        let theta2 = Math.PI * (i+1) / stacks;

        for (let j = 0; j < slices; j++) {
            let phi1 = 2 * Math.PI * j / slices;
            let phi2 = 2 * Math.PI * (j+1) / slices;

            // 4 corners of the patch
            let p1 = sphericalToCartesian(radius, theta1, phi1);
            let p2 = sphericalToCartesian(radius, theta2, phi1);
            let p3 = sphericalToCartesian(radius, theta2, phi2);
            let p4 = sphericalToCartesian(radius, theta1, phi2);

            // two triangles: (p1,p2,p3) and (p1,p3,p4)
            verts.push(...p1, ...p2, ...p3);
            verts.push(...p1, ...p3, ...p4);
        }
    }
    return verts;
}

function sphericalToCartesian(r, theta, phi) {
    return [
        r * Math.sin(theta) * Math.cos(phi),
        r * Math.cos(theta),
        r * Math.sin(theta) * Math.sin(phi),
    ];
}

function vertsToTriangles(verts) {
    let arr = [];

    for (let i = 0; i < verts.length/9; i++) {
        arr.push({
            v0: vec3(verts[i*9+0], verts[i*9+1], verts[i*9+2]),
            v1: vec3(verts[i*9+3], verts[i*9+4], verts[i*9+5]),
            v2: vec3(verts[i*9+6], verts[i*9+7], verts[i*9+8]),
            mat: createMaterial(vec3(0, 0, 0.5))
        });
    }

    return arr;
}

// quick helper to set a material color on a triangle array
function paint(tris, colorVec3, reflectivity=0, isPortal=false) {
    for (const t of tris) t.mat = createMaterial(colorVec3, reflectivity);
    return tris;
}

// Each face defined by 4 corners in CCW order as seen from outside,
// then split into two triangles: (0,1,2) and (0,2,3).
function quadToVerts(a, b, c, d) {
    return [
        ...a, ...b, ...c,   // tri 1
        ...a, ...c, ...d    // tri 2
    ];
}

// Corners (±H on each axis)
const nXnYnZ = [-H, -H, -H], nXnYpZ = [-H, -H,  H],
      nXpYnZ = [-H,  H, -H], nXpYpZ = [-H,  H,  H],
      pXnYnZ = [ H, -H, -H], pXnYpZ = [ H, -H,  H],
      pXpYnZ = [ H,  H, -H], pXpYpZ = [ H,  H,  H];

// +X face (x = +H), outward normal +X
const vertsPosX = quadToVerts(
  pXnYnZ, pXnYpZ, pXpYpZ, pXpYnZ
);
// -X face (x = -H), outward normal -X
const vertsNegX = quadToVerts(
  nXnYpZ, nXnYnZ, nXpYnZ, nXpYpZ
);

// +Y face (y = +H), outward normal +Y
const vertsPosY = quadToVerts(
  nXpYnZ, pXpYnZ, pXpYpZ, nXpYpZ
);
// -Y face (y = -H), outward normal -Y
const vertsNegY = quadToVerts(
  nXnYnZ, nXnYpZ, pXnYpZ, pXnYnZ
);

// +Z face (z = +H), outward normal +Z
const vertsPosZ = quadToVerts(
  nXnYpZ, nXpYpZ, pXpYpZ, pXnYpZ
);
// -Z face (z = -H), outward normal -Z
const vertsNegZ = quadToVerts(
  nXpYnZ, nXnYnZ, pXnYnZ, pXpYnZ
);

// Convert verts → triangles, then paint each face a different color
const trisPosX = paint(vertsToTriangles(vertsPosX), vec3(a, 0, 0), t); // red
const trisNegX = paint(vertsToTriangles(vertsNegX), vec3(0, a, 0), t); // green
const trisPosY = paint(vertsToTriangles(vertsPosY), vec3(0, 0, a), t); // blue
const trisNegY = paint(vertsToTriangles(vertsNegY), vec3(a, a, 0), t); // yellow
const trisPosZ = paint(vertsToTriangles(vertsPosZ), vec3(a, 0, a), t, true); // magenta
const trisNegZ = paint(vertsToTriangles(vertsNegZ), vec3(0, a, a), t); // cyan

// Final cube triangle list
const cubeTris = [
  ...trisPosX, 
  ...trisNegX,
  ...trisPosY, 
  ...trisNegY,
  ...trisPosZ, 
  ...trisNegZ
];

const triangles = [...cubeTris];

function vec3(x, y, z) {
    return {x: x, y: y, z: z};
}

function vec3add(v1, v2) {
    return vec3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

/** return `v1 - v2` */
function vec3sub(v1, v2) {
    return vec3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

function vec3dot(v1, v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

function vec3cross(v1, v2) {
    return vec3(
        v1.y * v2.z - v1.z * v2.y, // x component
        v1.z * v2.x - v1.x * v2.z, // y component
        v1.x * v2.y - v1.y * v2.x  // z component
    );
}

function vec3mul(v1, v2) {
    return vec3(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}

function vec3scale(k, v) {
    return vec3(v.x * k, v.y * k, v.z * k);
}

function normalize(v) {
    let mag = vec3mag(v);
    return vec3(v.x/mag, v.y/mag, v.z/mag);
}

function vec3mag(v) {
    return Math.hypot(v.x, v.y, v.z);
}

function vec3dist(v1, v2) {
    return vec3mag(vec3sub(v1, v2));
}

function rotateAroundAxis(v, axis, angle) {
    let cosA = Math.cos(angle);
    let sinA = Math.sin(angle);
    axis = normalize(axis);
    return vec3add(
        vec3add(
            vec3scale(cosA, v),
            vec3scale(sinA, vec3cross(axis, v))
        ),
        vec3scale((1 - cosA) * vec3dot(axis, v), axis)
    );
}

function rotateCamera(yaw, pitch) {
    // yaw: rotate around global/camera up
    camera.forward = rotateAroundAxis(camera.forward, camera.up, yaw);

    // pitch: rotate around camera right
    camera.right = normalize(vec3cross(camera.forward, camera.up));
    camera.forward = rotateAroundAxis(camera.forward, camera.right, pitch);

    // recompute orthonormal basis
    camera.right = normalize(vec3cross(camera.forward, camera.up));
    camera.up = normalize(vec3cross(camera.right, camera.forward));
}

rotateCamera(0, 0);

function makeRayFromCameraThroughPixel(x, y) {
    // idea: look into point in front of camera in 3D space (called the view plane). each pixel = spot on plane. origin = camera position, direction = camera pos to spot on view plane

    let u = ((x + 0.5) / WIDTH) * 2 - 1;
    let v = ((y + 0.5) / HEIGHT) * 2 - 1;

    // aspect ratio
    u *= WIDTH / HEIGHT;

    let planeHalfHeight = Math.tan(fov / 2);
    let planeHalfWidth = planeHalfHeight * WIDTH / HEIGHT;

    // finds point on view plane in world space
    let viewPlaneCenter = vec3add(camera.position, camera.forward);
    let pixelPoint = vec3add(viewPlaneCenter, vec3scale(u * planeHalfWidth, camera.right));
    pixelPoint = vec3sub(pixelPoint, vec3scale(v * planeHalfHeight, camera.up));

    // create ray, origin camera pos and dir to pixelPoint in world space calculated from view plane

    return {
        origin: camera.position,
        dir: normalize(vec3sub(pixelPoint, camera.position))
    };
}

// color: vec3, add other things later
function createMaterial(color, reflectivity=0, emissive = vec3(0, 0, 0)) {
    return {color: color, reflectivity, reflectivity, emissive};
}

// this would be a struct in c code
function createHitData(t, position, normal, obj_mat) {
    return {t: t, position: position, normal: normal, mat: obj_mat};
}

function intersectSphere(ray, sphere) {
    // ...
    let oc = vec3sub(ray.origin, sphere.position);
    let a = vec3dot(ray.dir, ray.dir); // should be 1 if normalized
    let b = 2 * vec3dot(oc, ray.dir);
    let c = vec3dot(oc, oc) - sphere.radius**2;

    let discriminant = b*b - 4*a*c;

    if (discriminant < 0) {
        return null; // object and ray cannot hit -- obj is behind ray/not in ray's path
    }

    let sqrtDisc = Math.sqrt(discriminant);
    let t1 = (-b - sqrtDisc) / (2*a);
    let t2 = (-b + sqrtDisc) / (2*a);

    let t; // distance to hit. figure out which one is the valid one. if neither are valid then that means the obj is behind the ray i think

    if (t1 > 0) t = t1;
    else if (t2 > 0) t = t2;
    else return null;

    let hitPoint = vec3add(ray.origin, vec3scale(t, ray.dir));
    let normal = normalize(vec3sub(hitPoint, sphere.position)); // wait this is so smart

    return createHitData(t, hitPoint, normal, sphere.mat);
}

function intersectTriangle(ray, triangle) {
    // O + tD = (1-u-v)V0 + uV1 + vV2

    let v0 = triangle.v0;
    let v1 = triangle.v1;
    let v2 = triangle.v2;

    let edge1 = vec3sub(v1, v0);
    let edge2 = vec3sub(v2, v0);

    let h = vec3cross(ray.dir, edge2);
    let a = vec3dot(edge1, h); // scalar triple product

    // u * (v x w) = 0 --> u parallel to v, w plane

    if (Math.abs(a) < 1e-6) {
        return null; // ray parallel to triangle
    }

    let f = 1.0 / a; // be careful in glsl with the integer division
    let s = vec3sub(ray.origin, v0);
    let u = f * vec3dot(s, h);

    if (u < 0 || u > 1) {
        return null; // outside triangle
    }

    let q = vec3cross(s, edge1);
    let v = f * vec3dot(ray.dir, q);

    if (v < 0 || (u + v) > 1) {
        return null;
    }

    let t = f * vec3dot(edge2, q);

    if (t > 1e-6) { // intersection in front of ray
        let hitPoint = vec3add(ray.origin, vec3scale(t, ray.dir));
        let normal = normalize(vec3cross(edge1, edge2));
        return createHitData(t, hitPoint, normal, triangle.mat);
    }
    else return null;
}

// finds two points where ray intersects with object, gives the dist closest to view point
function intersect(ray, obj) { // must be flexible for both spheres and rectangles
    if (obj.radius !== undefined) return intersectSphere(ray, obj);
    else if (obj.v0 !== undefined) return intersectTriangle(ray, obj);
    else {
        console.log(obj);
        throw new Error("not sphere or rect whats going on");
    }
}

// for (let i = 0; i < 10; i++) {
//     objs.push({
//         position: vec3(Math.random() * (20) - 10, Math.random() * (20) - 10, Math.random() * (20) - 10),
//         radius: 2,
//         mat: createMaterial(vec3(Math.random() * 0.75, 0.75 * Math.random(), 0.75 * Math.random()))
//     })
// }

// this is for finding the material the ray hits
function findClosestIntersection(ray) {
    let closestHit = null;
    let closestDist = Infinity;

    for (const obj of objs) {
        let hitData = intersect(ray, obj);
        if (hitData === null) continue; // in c i should do a isValid value in the hitData struct and check that

        let dist = hitData.t;

        if (dist > 0 && dist < closestDist) {
            closestHit = hitData;
            closestDist = dist;
        }
    }

    for (const triangle of triangles) {
        let hitData = intersect(ray, triangle);
        if (hitData === null) continue; // in c i should do a isValid value in the hitData struct and check that

        let dist = hitData.t;

        if (dist > 0 && dist < closestDist) {
            closestHit = hitData;
            closestDist = dist;
        }
    }

    // this is lowk so redundant but i write code that works, not code that looks good
    for (const light of lights) {
        let hitData = intersect(ray, light);
        if (hitData === null) continue; // in c i should do a isValid value in the hitData struct and check that

        let dist = hitData.t;

        if (dist > 0 && dist < closestDist) {
            closestHit = hitData;
            closestDist = dist;
        }
    }

    return closestHit; // could optimize this with space partioning but thats not necessary for this little mfs
}

function isInShadow(ray, light) {
    let ray_light_dist = vec3dist(ray.origin, light.position);
    for (const obj of objs) {
        let hitData = intersect(ray, obj);

        if (hitData == null) continue; // should be hitdata.isValid == false but that can be fixed later
        
        let dist = hitData.t;

        if (dist >= ray_light_dist) continue; // object behind light doesnt count
        else if (dist < ray_light_dist) return true; // object between light and ray, is in shadow, return true
    }
    
    return false;
}

function makeRay(origin, dir) {
    return {
        origin: origin,
        dir: dir
    }
}

function traceRay(ray, depth, allowPortal=false) {
    if (depth > MAX_DEPTH) return backgroundColor;
    if (depth <= 0) return vec3(0, 0, 0); // limit of reflection

    let hit = findClosestIntersection(ray);

    if (hit == null) return backgroundColor;

    let point = hit.position;
    let normal = hit.normal;
    let mat = hit.mat;

    if (mat.isPortal) {
        // Skip this hit: continue ray just past the wall
        let newOrigin = vec3add(hit.position, vec3scale(1e-4, ray.dir));
        return traceRay(makeRay(newOrigin, ray.dir), depth, false);
    }

    if (mat.emissive && (mat.emissive.x > 0 || mat.emissive.y > 0 || mat.emissive.z > 0)) {
        return mat.emissive;
    }

    // start with ambient light, then add contrib from every light source
    let finalColor = vec3mul(ambientLight, mat.color);

    for (const light of lights) {
        // add small episilon along normal to avoid ray hitting own object
        // ray from hit point to light
        let shadowRayOrigin = vec3add(point, vec3scale(1e-6, normal));
        let shadowRayDir = normalize(vec3sub(light.position, shadowRayOrigin));
        let shadowRay = makeRay(shadowRayOrigin, shadowRayDir);

        if (isInShadow(shadowRay, light)) continue; // object is in shadow, do not account in final lighting

        // angle based diffusion strength
        let diffuseStrength = Math.max(0, vec3dot(normal, normalize(vec3sub(light.position, point))));
        finalColor = vec3add(finalColor, vec3scale(diffuseStrength, vec3mul(light.color, mat.color)));
    }

    // reflection stuff here
    if (mat.reflectivity === undefined) return finalColor;
    if (mat.reflectivity > 0) {
        // reflection direction
        let reflectDir = vec3sub(ray.dir, vec3scale(2 * vec3dot(ray.dir, normal), normal));
        reflectDir = normalize(reflectDir);

        let reflectRay = makeRay(vec3add(point, vec3scale(1e-6, normal)), reflectDir);

        let reflectedColor = traceRay(reflectRay, depth - 1);

        finalColor = vec3add(vec3scale(1 - mat.reflectivity, finalColor), vec3scale(mat.reflectivity, reflectedColor));
    }

    return finalColor;
}

for (let y = 0; y < canvas.height; y++) {
    for (let x = 0; x < canvas.width; x++) {
        let ray = makeRayFromCameraThroughPixel(x, y);
        let color = traceRay(ray, MAX_DEPTH, true);

        let r = color.x * 255;
        let g = color.y * 255;
        let b = color.z * 255;
        ctx.fillStyle = `rgb(${r}, ${g}, ${b})`;
        ctx.fillRect(x, y, 1, 1);
    }
}
